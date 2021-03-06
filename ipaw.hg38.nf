/*
vim: syntax=groovy
-*- mode: groovy;-*-

==============================
IPAW: HiRIEF II varDB pipeline
==============================
@Authors
Jorrit Boekel @glormph
Yafeng Zhu @yafeng

https://github.com/lehtiolab/proteogenomics-analysis-workflow

- make sure conda install all excutive programs or scripts first
- create a bin folder to store all scripts
*/

nf_required_version = '0.26.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}


/* SET DEFAULT PARAMS */
params.isobaric = false
params.bamfiles = false
params.outdir = 'result'

params.novheaders = '^PGOHUM;^lnc;^decoy_PGOHUM;^decoy_lnc' 
params.varheaders = '^COSMIC;^CanProVar;^decoy_COSMIC;^decoy_CanProVar'
params.saavheader = false
params.noclassfdr = false

mods = file(params.mods)
knownproteins = file(params.knownproteins)
blastdb = file(params.blastdb)
gtffile = file(params.gtf)
snpfa = file(params.snpfa)

genomefa = file(params.genome)
msgf_db = file(params.tdb)
msgf_jar = file(params.msgfjar)
converter = file(params.mzid2tsvConverter)

activationtype = [HCD:'High-energy collision-induced dissociation',CID:'Collision-induced dissociation',ETD:'Electron transfer dissociation',no:''][params.activationFilter]
massshifts = [tmt:0.0013, itraq:0.00125, false:0]
plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : false
massshift = massshifts[plextype]
msgfprotocol = params.isobaric ? [tmt:4, itraq:2][plextype] : 0

params.PrecursorMassTolerance = '10ppm'
params.FragmentMethodID = 3 // 3 for HCD spectra, 1 for CID, 2 for ETD


params.tasks = 1
params.thread = 1
params.qval = 0.05 
params.MS2error = 0.02 // unit: Da
//////////////////
// FOR NON-VARDB DATABASES
// E.g. 6FT, 3FT or WXS
novheaders = params.novheaders == true ? false : params.novheaders
varheaders = params.varheaders == true ? false : params.varheaders
///////////////////

println "PrecursorMassTolerance used in MSGF+:${params.PrecursorMassTolerance}"
println "q-value cutoff:${params.qval}"
println "MS2 error tolerance used in SpectrumAI:${params.MS2error}"
println "IsobaricAnalyzer filtered by activation:${params.activationFilter}"
println "InstrumentID used in MSGF+: ${params.inst}"
println "FragmentMethodID used in MSGF+: ${params.FragmentMethodID}"
println "ProtocolID used in MSGF+: $msgfprotocol"


/* PIPELINE START */

// Either feed an mzmldef file (tab separated lines with filepath\tsetname), or /path/to/\*.mzML
if (!params.mzmldef) {
Channel
  .fromPath(params.mzmls)
  .map { it -> [it, 'NA'] }
  .set { mzml_in }
} else {
Channel
  .from(file("${params.mzmldef}").readLines())
  .map { it -> it.tokenize('\t') }
  .set { mzml_in }
}


mzml_in
  .tap { sets }
  .map { it -> [ it[1], (it[0]=~/.*\/(.+)\..*$/)[0][1], file(it[0])]} // create, set, samplename and file
  .tap{ mzmlfiles; mzml_msgf ; mzml_isobaric }
  .count()
  .set{ amount_mzml }

sets
  .map{ it -> it[1] }
  .unique()
  .tap { sets_for_emtpybam; sets_for_denoms }
  .collect()
  .subscribe { println "Detected setnames: ${it.join(', ')}" }


// Get denominators if isobaric experiment
// passed in form --denoms 'set1:126:128N set2:131 set4:129N:130C:131'
if (params.isobaric && params.denoms) {
  setdenoms = [:]
  params.denoms.tokenize(' ').each{ it -> x=it.tokenize(':'); setdenoms.put(x[0], x[1..-1])}
  set_denoms = Channel.value(setdenoms)  
}else if (params.isobaric) {
  setdenoms = [:]
  sets_for_denoms.reduce(setdenoms){ a, b -> a.put(b, ['_126']); return a }.set{ set_denoms }
}

process makeProtSeq {

  input:
  file knownproteins

  output:
  file('mslookup_db.sqlite') into protseqdb

  """
  msslookup protspace -i $knownproteins --minlen 8
  """
}

process makeTrypSeq {

  input:
  file knownproteins

  output:
  file('mslookup_db.sqlite') into trypseqdb

  """
  msslookup seqspace -i $knownproteins --insourcefrag
  """
}


process IsobaricQuant {

  when: params.isobaric

  input:
  set val(setname), val(sample), file(infile) from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  """
  IsobaricAnalyzer  -type $params.isobaric -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
  """
}

isobaricamount = params.isobaric ? amount_mzml.value : 1

isobaricxml
  .ifEmpty(['NA', 'NA'])
  .buffer(size: isobaricamount)
  .flatMap { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> it[1] }
  .collect()
  .set { sorted_isoxml }

mzmlfiles
  .tap { groupset_mzmls }
  .buffer(size: amount_mzml.value)
  .map { it.sort( {a, b -> a[1] <=> b[1]}) }  //to sort by the first element of a tuple in descending order
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }] }
  .set{ mzmlfiles_all }

process createSpectraLookup {

  input:
  file(isobxmls) from sorted_isoxml 
  set val(setnames), file(mzmlfiles) from mzmlfiles_all
  
  output:
  file('mslookup_db.sqlite') into spec_lookup

  script:
  if(params.isobaric)
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  msslookup isoquant --dbfile mslookup_db.sqlite -i ${isobxmls.join(' ')} --spectra ${mzmlfiles.join(' ')}
  """
  else
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}


process msgfPlus {

  // Latest version has problems when converting to TSV, possible too long identifiers 
  // So we use an older version 2016.10.26
  // LATEST TESTED: container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'  

  input:
  set val(setname), val(sample), file(x) from mzml_msgf

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  
  """
  fs=`du -Lk $msgf_db|cut -f1`
  java -Xmx\$((\$fs*8/1024))M -jar $msgf_jar -thread ${params.thread} -tasks ${params.tasks} -d $msgf_db -s $x -o "${sample}.mzid" -mod $mods -tda 0 -t ${params.PrecursorMassTolerance} -ti -1,2 -m ${params.FragmentMethodID} -inst ${params.inst} -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 7 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1
  """
}

mzids
  .tap { mzids_2tsv }
  .groupTuple()
  .set { mzids_2pin }

process ConvertmzidTotsv {

  input:
  set val(setname), val(sample), file(x) from mzids_2tsv
  
  output:
  set val(setname), file("${sample}.mzid"), file('out.mzid.tsv') into mzidtsvs 

  """
  mono $converter -mzid:$x -tsv:temp.tsv
  cut_extra_aa.py temp.tsv out.mzid.tsv
  rm -f temp.tsv
  """
}

process percolator {

  input:
  set val(setname), val(samples), file('mzid?') from mzids_2pin

  output:
  set val(setname), file('perco.xml') into percolated

  """
  echo $samples
  mkdir mzids
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/mzid\$count mzids/\${sam}.mzid; echo mzids/\${sam}.mzid >> metafile; ((count++));done
  msgf2pin -o percoin.xml -e trypsin -P "decoy_" metafile
  percolator -j percoin.xml -X perco.xml -N 500000 --decoy-xml-output -y
  """
}


percolated
  .tap { var_percolated }
  .set { nov_percolated }


process getVariantPercolator {

  when: varheaders

  input:
  set val(setname), file(x) from var_percolated

  output:
  set val(setname), val('variant'), file("${x}_h0.xml") into var_perco
  """
  msspercolator splitprotein -i $x --protheaders \'$varheaders\'
  """
}


process getNovelPercolator {

  when: novheaders

  input:
  set val(setname), file(x) from nov_percolated

  output:
  set val(setname), val('novel'), file("${x}_h0.xml") into nov_perco
  """
  msspercolator splitprotein -i $x --protheaders \'$novheaders\'
  """
}


nov_perco
  .concat(var_perco)
  .set { splitperco }


process filterPercolator {

  input:
  set val(setname), val(peptype), file(perco) from splitperco
  file 'trypseqdb' from trypseqdb
  file 'protseqdb' from protseqdb
  file knownproteins

  output:
  set val(setname), val(peptype), file('filtprot') into filtered_perco
  
  script:
  if (params.noclassfdr)
  """
  mv $perco filtprot
  """
  else
  """
  msspercolator filterseq -i $perco -o filtseq --dbfile trypseqdb --insourcefrag 2 --deamidate 
  msspercolator filterprot -i filtseq -o filtprot --fasta $knownproteins --dbfile protseqdb --minlen 8 --deamidate --enforce-tryptic
  """
}

nov_filtered_perco = Channel.create()
var_filtered_perco = Channel.create()
filtered_perco
  .choice( var_filtered_perco, nov_filtered_perco) { it -> it[1] == 'variant' ? 0 : 1 }

mzidtsvs
  .groupTuple()
  .tap { variantmzidtsv }
  .join(nov_filtered_perco)
  .set { nov_mzperco }

variantmzidtsv
  .join(var_filtered_perco)
  .concat(nov_mzperco)
  .set { allmzperco }

process svmToTSV {

  input:
  set val(setname), file('mzident?'), file('mzidtsv?'), val(peptype), file(perco) from allmzperco 

  output:
  set val(setname), val(peptype), file('mzidperco') into mzidtsv_perco

  script:
  """
#!/usr/bin/env python3
from glob import glob
mzidtsvfns = sorted(glob('mzidtsv*'))
mzidfns = sorted(glob('mzident*'))
from app.readers import pycolator, xml, tsv, mzidplus
import os
ns = xml.get_namespace_from_top('$perco', None) 
psms = {p.attrib['{%s}psm_id' % ns['xmlns']]: p for p in pycolator.generate_psms('$perco', ns)}
decoys = {True: 0, False: 0}
for psm in sorted([(pid, float(p.find('{%s}svm_score' % ns['xmlns']).text), p) for pid, p in psms.items()], reverse=True, key=lambda x:x[1]):
    pdecoy = psm[2].attrib['{%s}decoy' % ns['xmlns']] == 'true'
    decoys[pdecoy] += 1
    try:
        psms[psm[0]] = {'decoy': pdecoy, 'svm': psm[1], 'qval': decoys[True]/decoys[False]}  # T-TDC
    except ZeroDivisionError:
        psms[psm[0]] = {'decoy': pdecoy, 'svm': psm[1], 'qval': 1.0}  # T-TDC  
decoys = {'true': 0, 'false': 0}
for svm, pep in sorted([(float(x.find('{%s}svm_score' % ns['xmlns']).text), x) for x in pycolator.generate_peptides('$perco', ns)], reverse=True, key=lambda x:x[0]):
    decoys[pep.attrib['{%s}decoy' % ns['xmlns']]] += 1
    try:
        [psms[pid.text].update({'pepqval': decoys['true']/decoys['false']}) for pid in pep.find('{%s}psm_ids' % ns['xmlns'])]
    except ZeroDivisionError:
        [psms[pid.text].update({'pepqval': 1.0}) for pid in pep.find('{%s}psm_ids' % ns['xmlns'])]   
oldheader = tsv.get_tsv_header(mzidtsvfns[0])
header = oldheader + ['percolator svm-score', 'PSM q-value', 'peptide q-value']
with open('mzidperco', 'w') as fp:
    fp.write('\\t'.join(header))
    for fnix, mzidfn in enumerate(mzidfns):
        mzns = mzidplus.get_mzid_namespace(mzidfn)
        siis = (sii for sir in mzidplus.mzid_spec_result_generator(mzidfn, mzns) for sii in sir.findall('{%s}SpectrumIdentificationItem' % mzns['xmlns']))
        for specidi, psm in zip(siis, tsv.generate_tsv_psms(mzidtsvfns[fnix], oldheader)):
            # percolator psm ID is: samplename_SII_scannr_rank_scannr_charge_rank
            print(specidi)
            print(psm)
            scan, rank = specidi.attrib['id'].replace('SII_', '').split('_')
            outpsm = {k: v for k,v in psm.items()}
            spfile = os.path.splitext(psm['#SpecFile'])[0]
            try:
                percopsm = psms['{fn}_SII_{sc}_{rk}_{sc}_{ch}_{rk}'.format(fn=spfile, sc=scan, rk=rank, ch=psm['Charge'])]
            except KeyError:
                continue
            if percopsm['decoy']:
                continue
            fp.write('\\n')
            outpsm.update({'percolator svm-score': percopsm['svm'], 'PSM q-value': percopsm['qval'], 'peptide q-value': percopsm['pepqval']})
            fp.write('\\t'.join([str(outpsm[k]) for k in header]))
  """
}

mzidtsv_perco
  .combine(spec_lookup)
  .set { prepsm }

process createPSMTables {

  input:
  set val(setname), val(peptype), file('psms'), file('lookup') from prepsm

  output:
  set val(setname), val(peptype), file("${setname}_${peptype}_psmtable.txt") into psmtable

  script:
  if(params.isobaric)
  """
  msspsmtable conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl ${params.qval} --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl ${params.qval} --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup
  msspsmtable specdata -i filtpep --dbfile psmlookup -o prepsms.txt
  msspsmtable quant -i prepsms.txt -o ${setname}_${peptype}_psmtable.txt --dbfile psmlookup --isobaric
  sed 's/\\#SpecFile/SpectraFile/' -i ${setname}_${peptype}_psmtable.txt
  """
  else
  """
  msspsmtable conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl ${params.qval} --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl ${params.qval} --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup
  msspsmtable specdata -i filtpep --dbfile psmlookup -o ${setname}_${peptype}_psmtable.txt
  sed 's/\\#SpecFile/SpectraFile/' -i ${setname}_${peptype}_psmtable.txt
  """
}


variantpsms = Channel.create()
novelpsms = Channel.create()
psmtable
  .tap { setmergepsmtables; peppsms }
  .choice( variantpsms, novelpsms ) { it -> it[1] == 'variant' ? 0 : 1 }

 
setmergepsmtables
  .groupTuple(by: 1)
  .set { psmmerge_in }


process mergeSetPSMtable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setnames), val(peptype), file(psmtables) from psmmerge_in

  output:
  file "${peptype}_psmtable.txt" into produced_psmtables

  """
  head -n1 ${psmtables[0]} > ${peptype}_psmtable.txt
  for fn in ${psmtables.join(' ')}; do tail -n+2 \$fn >> ${peptype}_psmtable.txt; done
  """
}


process prePeptideTable {
  
  input:
  set val(setname), val(peptype), file('psms?') from peppsms 

  output:
  set val(setname), val(peptype), file('peptidetable.txt') into peptable
 
  script:
  """
  msspsmtable merge -o psms.txt -i psms*
  msspeptable psm2pep -i psms.txt -o preisoquant --scorecolpattern svm --spectracol 1 ${params.isobaric ? '--isobquantcolpattern plex' : ''}
  awk -F '\\t' 'BEGIN {OFS = FS} {print \$12,\$13,\$3,\$7,\$8,\$9,\$11,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22}' preisoquant > preordered
  ${params.isobaric ? "msspsmtable isoratio -i psms.txt -o peptidetable.txt --targettable preordered --isobquantcolpattern plex --minint 0.1 --denompatterns ${set_denoms.value[setname].join(' ')} --protcol 11" : 'mv preordered peptidetable.txt'}    

  """
}

novelpsms
  .into{novelpsmsFastaBedGFF; novelpsms_specai}

novelprepep = Channel.create()
presai_peptable = Channel.create()
peptable
  .choice( presai_peptable, novelprepep ) { it -> it[1] == 'variant' ? 0 : 1 }
novelprepep
  .join(novelpsmsFastaBedGFF)
  .set { novelFaBdGfPep }

process createFastaBedGFF {
 
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "${setname}_novel_peptides.gff3" ? "${setname}_novel_peptides.gff3" : null}
 
  input:
  set val(setname), val(peptype), file(peptides) , val(psmtype), file(psms) from novelFaBdGfPep  

  output:
  set val(setname), file('novel_peptides.fa') into novelfasta
  set val(setname), file('novel_peptides.bed') into novelbed
  set val(setname), file("${setname}_novel_peptides.gff3") into novelGFF3
  set val(setname), file('novel_peptides.tab.txt') into novelpep
  set val(setname), file('novpep_perco_quant.txt') into novelpep_percoquant
 
  """
  map_novelpeptide2genome.py --input $psms --gtf $gtffile --fastadb $msgf_db --tab_out novel_peptides.tab.txt --fasta_out novel_peptides.fa --gff3_out ${setname}_novel_peptides.gff3 --bed_out novel_peptides.bed
  sort -k 1b,1 <(tail -n+2 $peptides) |cut -f 1,14-500 > peptable_sorted
  sort -k 2b,2 <(tail -n+2 novel_peptides.tab.txt) > novpep_sorted
  paste <(cut -f 2 novpep_sorted) <(cut -f1,3-500 novpep_sorted) > novpep_pepcols
  join novpep_pepcols peptable_sorted -j 1 -a1 -o auto -e 'NA' -t \$'\\t' > novpep_pqjoin
  paste <(cut -f 2 novpep_pqjoin) <(cut -f1,3-500 novpep_pqjoin) > novpep_joined_pepcols
  paste <(head -n1 novel_peptides.tab.txt)  <(cut -f 14-500 $peptides |head -n1) > header
  cat header novpep_joined_pepcols > novpep_perco_quant.txt
  """
}

novelpep
  .into {blastnovelpep; blatnovelpep; annonovelpep; snpnovelpep}
novelfasta
  .into {blastnovelfasta; blatnovelfasta}

process BlastPNovel {

  input:
  set val(setname), file(novelfasta) from blastnovelfasta
  file blastdb

  output:
  set val(setname), file('blastp_out.txt') into novelblast
  
  """
  makeblastdb -in $blastdb -dbtype prot
  blastp -db $blastdb -query $novelfasta -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 4 -max_target_seqs 1 -evalue 1000 -out blastp_out.txt
  """
}

novelpsms_specai
  .map { it -> [it[0], it[2]] }
  .join(blastnovelpep)
  .join(novelblast)
  .set { novelblastout }

process ParseBlastpOut {

 input:
 set val(setname), file(psms), file(novelpep), file(novelblast) from novelblastout
 file blastdb

 output:
 set val(setname), file('peptable_blastp.txt') into peptable_blastp
 set val(setname), file('single_mismatch_novpeps.txt') into novpeps_singlemis

 """
 parse_BLASTP_out.py --input $novelpep --blastp_result $novelblast --fasta $blastdb --output peptable_blastp.txt
 extract_1mismatch_novpsm.py peptable_blastp.txt $psms single_mismatch_novpeps.txt
 """

}

groupset_mzmls
  .map { it -> [it[0], it[1], it[2]] } // strip fraction, strip
  .groupTuple()
  .tap { var_specaimzmls }
  .join(novpeps_singlemis)
  .set { grouped_saavnov_mzml_peps }
  

process ValidateSingleMismatchNovpeps {
 
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "precursorError.histogram.plot.pdf" ? "${setname}_novel_precursorError_plot.pdf" : it }
  
  input:
  set val(setname), val(samples), file(mzmls), file(peps) from grouped_saavnov_mzml_peps

  output:
  set val(setname), file("${setname}_novel_saav_specai.txt") into singlemis_specai


  """
  mkdir mzmls
  for fn in $mzmls; do ln -s `pwd`/\$fn mzmls/; done
  Rscript ../../../SpectrumAI/SpectrumAI.R mzmls $peps ${setname}_novel_saav_specai.txt ${params.MS2error} || cp $peps singlemis_specai.txt
  """
}

singlemis_specai
  .join(peptable_blastp)
  .set { nov_specaiparse }

process novpepSpecAIOutParse {

  input:
  set val(setname), file(x), file('peptide_table.txt') from nov_specaiparse
  
  output:
  set val(setname), file('novpep_specai.txt') into novpep_singlemisspecai

  """
  parse_spectrumAI_out.py --spectrumAI_out $x --input peptide_table.txt --output novpep_sa
  cut -f 1,8-19 novpep_sa > novpep_specai.txt
  """
}


process BLATNovel {

  input:
  set setname, file(novelfasta) from blatnovelfasta
  file genomefa

  output:
  set val(setname), file('blat_out.pslx') into novelblat

  """
  blat $genomefa $novelfasta -t=dnax -q=prot -tileSize=5 -minIdentity=99 -out=pslx blat_out.pslx

  """
}

novelblat
  .join(blatnovelpep)
  .set { novblatparse }

process parseBLATout {

 input:
 set val(setname), file(novelblat), file(novelpep) from novblatparse

 output:
 set val(setname), file('peptable_blat.txt') into peptable_blat

 """
 parse_BLAT_out.py $novelblat $novelpep peptable_blat.txt

 """
}

process labelnsSNP {
  
  input:
  set val(setname), file(peptable) from snpnovelpep
  file snpfa

  output:
  set val(setname), file('nssnp.txt') into ns_snp_out

  """
  label_nsSNP_pep.py --input $peptable --nsSNPdb $snpfa --output nssnp.txt
  """
}

novelGFF3
  .into { novelGFF3_phast; novelGFF3_phylo; novelGFF3_bams }

process phastcons {
 
  input:
  set val(setname), file(novelgff) from novelGFF3_phast
  output:
  set val(setname), file ('phastcons.txt') into phastcons_out

  """
  calculate_phastcons.py $novelgff ../../../hg38.phastCons100way.bw phastcons.txt
  """
}

process phyloCSF {

  input:
  set val(setname), file(novelgff) from novelGFF3_phylo

  output:
  set val(setname), file('phylocsf.txt') into phylocsf_out

  """
  calculate_phylocsf.py $novelgff ../../../PhyloCSFtracks_hg38 phylocsf.txt
  """

}


if (params.bamfiles) {
  bamFiles = Channel
    .fromPath(params.bamfiles)
    .map { fn -> [ fn, fn + '.bai' ] }
    .collect()
} else {
  bamFiles = Channel.empty()
}


process scanBams {

  when: params.bamfiles

  input:
  set val(setname), file(gff) from novelGFF3_bams
  file bams from bamFiles
  
  output:
  set val(setname), file('scannedbams.txt') into scannedbams

  """
  ls *.bam > bamfiles.txt
  scan_bams.py  --gff_input $gff --bam_files bamfiles.txt --output scannedbams.txt
  """
}


process annovar {
 
  input:
  set val(setname), file(novelbed) from novelbed

  output:
  set val(setname), file('novpep_annovar.variant_function') into annovar_out

  """
  annotate_variation.pl -out novpep_annovar -build hg38 $novelbed ../../../humandb_v2/
  """

}

annovar_out
  .join(annonovelpep)
  .set { parseanno }

process parseAnnovarOut {
  
  input:
  set val(setname), file(anno), file(novelpep) from parseanno

  output:
  set val(setname), file('parsed_annovar.txt') into annovar_parsed

  """
  parse_annovar_out.py --input $novelpep --output parsed_annovar.txt --annovar_out $anno 
  """
}

if (params.bamfiles){
  scannedbams
    .set { bamsOrEmpty }
}
else {
  sets_for_emtpybam
    .map { it -> [it, 'No bams'] }
    .set { bamsOrEmpty }
}

ns_snp_out
  .join(novpep_singlemisspecai)
  .join(peptable_blat)
  .join(annovar_parsed)
  .join(phastcons_out)
  .join(phylocsf_out)
  .join(novelpep_percoquant)
  .join(bamsOrEmpty)
  .set { combined_novel }


process combineResults{

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setname), file(a), file(b), file(c), file(d), file(e), file(f), file(g), file(h) from combined_novel
  
  output:
  set val('nov'), val(setname), file("${setname}_novel_peptides.txt") into novpeps_finished
  
  script:
  if (!params.bamfiles)
  """
  for fn in $a $b $c $d $e $f $g; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  join joined5 $g -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined6
  grep '^Peptide' joined6 > ${setname}_novel_peptides.txt
  grep -v '^Peptide' joined6 >> ${setname}_novel_peptides.txt
  """

  else
  """
  for fn in $a $b $c $d $e $f $g $h; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  join joined5 $g -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined6
  join joined6 $h -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined7
  grep '^Peptide' joined7 > ${setname}_novel_peptides.txt
  grep -v '^Peptide' joined7 >> ${setname}_novel_peptides.txt
  """
}

process prepSpectrumAI {
  
  input:
  set val(setname), val(peptype), file(psms) from variantpsms
  
  output:
  set val(setname), file('specai_in.txt') into var_specai_input
  
  script:
  if (params.saavheader)
  """
  cat <(head -n1 $psms) <(grep $params.saavheader $psms) > saavpsms
  label_sub_pos.py --input_psm saavpsms --output specai_in.txt ${params.splitchar ? "--splitchar ${params.splitchar}" : ''}
  """
  else
  """
  label_sub_pos.py --input_psm $psms --output specai_in.txt
  """
}


var_specaimzmls
  .join(var_specai_input)
  .set { var_specai_inmzml }

process SpectrumAI { 

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setname), val(samples), file(mzmls), file(specai_in) from var_specai_inmzml

  output:
  set val(setname), file("${setname}_variant_specairesult.txt") into specai

  """
  mkdir mzmls
  for fn in $mzmls; do ln -s `pwd`/\$fn mzmls/; done
  ls mzmls
  Rscript ../../../SpectrumAI/SpectrumAI.R mzmls $specai_in ${setname}_variant_specairesult.txt ${params.MS2error}
  """
}

specai
  .join(presai_peptable)
  .set { specai_peptable }

process generateVariantPepTable {
 
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setname), file(x), val(peptype), file(peptides) from specai_peptable
  
  output:
  set val('var'), val(setname), file("${setname}_variant_peptides.txt") into varpeps_finished

  """
  ${params.saavheader ? "cat <(head -n1 ${peptides}) <(grep ${params.saavheader} ${peptides}) > saavpeps" : "mv ${peptides} saavpeps" }
  parse_spectrumAI_out.py --spectrumAI_out $x --input saavpeps --output setsaavs
  ${params.saavheader ? "cat setsaavs <(grep -v ${params.saavheader} ${peptides} | sed \$'s/\$/\tNA/') > ${setname}_variant_peptides.txt" : "mv setsaavs ${setname}_variant_peptides.txt"}
 
  # Remove PSM-table specific stuff (RT, precursor, etc etc) from variant PEPTIDE table
  cut -f 1,2,14-5000 ${setname}_variant_peptides.txt > pepsfix
  mv pepsfix ${setname}_variant_peptides.txt
  """
}

novpeps_finished
  .concat(varpeps_finished) 
  .groupTuple()
  .set { setmerge_peps }

accession_keymap = ['var': 'Peptide sequence', 'nov': 'Mod.peptide']
acc_removemap = ['nov': 'Peptide', 'var': 'Mod.peptide']

process mergeSetPeptidetable {
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(peptype), val(setnames), file('peps?') from setmerge_peps

  output:
  file "${peptype}_peptidetable.txt" into produced_peptables

  """
  # build non-changing fields (seq based fields) table:
  fixfields=`head -n1 peps1 |tr -s '\\t' '\\n' | egrep -vn '(Setname|Spectrum|q-val|plex|${acc_removemap[peptype]})' | cut -f 1 -d ':'`
  fixfields=`echo \$fixfields | sed 's/ /,/g'`
  head -n1 peps1 | cut -f `echo \$fixfields` > fixheader
  count=1; for setn in ${setnames.join(' ')} ; do
    cut -f `echo \$fixfields` peps\$count | tail -n+2 >> fixpeps
    ((count++))
  done
  if [ ${peptype} == 'nov' ]
  then 
     cat fixheader <(sort -u -k1b,1 fixpeps) > temp
     group_novpepToLoci.py  --input temp --output temp.loci --distance 10kb
     head -n1 temp.loci > fixheader
     tail -n+2 temp.loci > fixpeps
  fi
  sort -u -k1b,1 fixpeps > temp
  mv temp fixpeps

  ## Build changing fields table
  touch peptable
  count=1; for setn in ${setnames.join(' ')}; do
    varfields=`head -n1 peps\$count |tr -s '\\t' '\\n' | egrep -n '(${accession_keymap[peptype]}|Spectrum|q-val|plex)' | cut -f 1 -d ':'`
    varfields=`echo \$varfields| sed 's/ /,/g'`
    # first add to header, cut from f2 to remove join-key pep seq field
    head -n1 peps\$count | cut -f `echo \$varfields` | cut -f 2-5000| sed "s/^\\(\\w\\)/\${setn}_\\1/;s/\\(\\s\\)/\\1\${setn}_/g" > varhead
    paste fixheader varhead > newheader && mv newheader fixheader
    # then join the values
    tail -n+2 peps\$count | cut -f `echo \$varfields` | sort -k1b,1 > sortpep; join peptable sortpep -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined
    mv joined peptable
    ((count++))
  done
  join fixpeps peptable -a1 -a2 -o auto -e 'NA' -t \$'\\t' > fixvarpeps
  cat fixheader fixvarpeps > ${peptype}_peptidetable.txt 
  """
}
produced_psmtables
  .concat(produced_peptables)
  .subscribe { println "Pipeline output ready: ${it}" }
