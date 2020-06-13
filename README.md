# pan-cancer-proteogenomics-analysis

This is the workflow of  the "Proteogenomics analysis of non-coding region encoded peptides in normal tissues and five cancer types" article. This workflow is modified from a previously published 
[Integrated proteogenomics analysis workflow](https://github.com/lehtiolab/proteogenomics-analysis-workflow) (IPAW)


# 1. LC-MS/MS Raw Files 

Prepare LC-MS/MS Raw Files. LC-MS/MS raw files of 40 normal samples from 31 healthy tissues and 926 cancer samples from five cancer types were obtained from [National Cancer Institute Clinical Proteomic Tumor Analysis Consortium ](<https://cptac-data-portal.georgetown.edu/cptac/public>) and [PRoteomics IDEntifcations](<https://www.ebi.ac.uk/pride/>) database

# 2. Data conversion

The IPAW requires the mzML files format as input files. Therefore, Convert Raw files to mzML files by [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser).

convert thr raw files using one core:
 (raw files are listed in rawfilelist.txt, one line for each file)
```bash
while read f;do mono ThermoRawFileParser.exe -i="$f" -o=/hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MS_Data/LiverCancer/ -f=1 -m=1; done < rawfilelist.txt
```
`-o` :the converted  output file path
` rawfilelist.txt` : contains the raw data absolute path.

to run the converting in multiple cores in parallel:

```
cat rawfilelist.txt |parallel -j 24 mono ThermoRawFileParser.exe -i={} -o=/data/home/yz332/MS_data/PXD000561_mzML/ -f=1 -m=1

```
`-o` :the converted  output file path
`-j` : the numbers of cores to run in parallel

CPTAC usually shares the MS raw data in mzML.gz format. To extract mzML files in specified folder, use:
```
for f in *Proteome*/*_mzML/*.gz; do 
STEM=$(basename "${f}" .gz)
gunzip -c "${f}" > /path/"${STEM}"
done
```

`/path/` : absolute file path for the extracted mzML files.

# 3. Database construction

**The healthy tissues**

To search the proteomics data of the healthy tissues, we used a core database which contains the human protein database of **Ensembl 92, CanProVar 2.0 variant peptides and peptide sequences from three frame translation of annotated pseudogenes and lncRNAs**.Pseudogenes were downloaded from GENCODE v28 including both annotated and predicted12. LncRNAs were downloaded from LNCpedia 4.1.

1. The core database : `VarDB2+Ensembl92.noCOSMIC.fa`

2. Build the decoy database: the decoy peptide was produced by reversing protein sequences in the target database.

   ```bash
   decoy.py VarDB2+Ensembl92.noCOSMIC.fa VarDB2+Ensembl92.noCOSMIC.revCat.fa
   # create index for the decoy databas
   java -Xmx10000M -cp /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MSGFplus/MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d VarDB2+Ensembl92.noCOSMIC.revCat.fa -tda 0
   ```

**The tumor tissues**
For each cancer type, **genomics data detected mutations** were downloaded from the [Cancer  Genomic Data Server](http://www.cbioportal.org/datasets) (CGDS) and then converted to mutant peptide sequences, which were added to the core database.
For example，lung cancer database construction

1. download the file [data_mutations_extended.txt](http://www.cbioportal.org/datasets)

2. Build a variation database

   ```bash
   Python3 mutations_to_proteindb.py lungcancer_path1/data_mutations_extended.txt Ensembl92+75.cds.all.fa lung_path1/mutproteins.fa
   Python3 mutations_to_proteindb.py lungcancer_path2/data_mutations_extended.txt Ensembl92+75.cds.all.fa(ref) lung_path2/mutproteins.v2.fa
   ...
   Python3 digest_mutant_protein.py Ensembl92+75.cds.all.fa lung_mutpeptides.txt lungcancer_path1/mutproteins.fa lungcancer_path2/mutproteins.fa lungcancer_path3/mutproteins.fa lungcancer_path4/mutproteins.fa ...
   ```

3. Add to the core database.

   ```bash
   cat /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/proteogenomics-analysis-workflow/VarDB2.0/VarDB2+Ensembl92.noCOSMIC.fa lung_mutpeptides.txt > VarDB2+lungMut+Ensembl92.fa
   ```

4. Build the decoy database
   The decoy peptide was produced by reversing protein sequences in the target database.

   ```bash
   Python3 decoy.py VarDB2+lungMut+Ensembl92.fa VarDB2+lungMut+Ensembl92.revCat.fa
   # create index for the decoy database
   java -Xmx10000M -cp MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d VarDB2+lungMut+Ensembl92.revCat.fa -tda 0
   ```

# 4. Analyse your mzML files

Customize the input files and parameter:
1. A tab deliminated text file with mzmlfilepath(absolute path) and setname 
   `--mzmldef  adrenal_gland_filelist_setnames.txt`

if you find out each set contains the same number of MS data file, you can quickly create filelist_setnames in this way:

for example, if there are in total 6 sets(experiments), each set contains 12 fractions/files:
```
ls -1 $PWD/*.mzML >filelist.txt
for i in {1..6};do printf "set$i"'%0.s\n' {1..12};done >setnames.txt
paste filelist.txt setnames.txt > filelist_setnames.txt
```

2. Quantitative method
   `--mods labelfree_Mods.txt`
    *     itraq4plex_Mods.txt  
    *     itraq8plex_Mods.txt  
    *     labelfree_Mods.txt  
    *     tmt_Mods.txt
3. Instrument 
   ` --inst 3`
    *    0: Low-res LCQ/LTQ (Default)
    *    1: Orbitrap/FTICR/Lumos
    *    2: TOF
    *    3: Q-Exactive

More custom parameters see [here](https://github.com/lehtiolab/proteogenomics-analysis-workflow)

```bash
#!/bin/bash
export NXF_WORK="/ldfssz1/ST_CANCER/CGR/USER/maleyao/IPAW_work/Normal/adrenal_gland/work"
path="/hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/proteogenomics-analysis-workflow/"
path2="/ldfssz1/ST_CANCER/CGR/USER/maleyao/mutationDB"

nextflow run ipaw.hg38.MSGFplus2018.v6.nf --msgfjar /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MSGFplus/MSGFPlus.jar \
--tdb $path2/VarDB2+Ensembl92.noCOSMIC.revCat.fa \ ## target decoy combined databases
--gtf $path/VarDB2.0/VarDB2.gtf \
--inst 3 \
--mods labelfree_Mods.txt \ #Modification file for MSGF+. Default file is for TMT labelled samples
--mzmldef adrenal_gland_filelist_setnames.txt # a tab deliminated text file with mzmlfilepath(absolute path) and setname 
--knownproteins $path/Ensembl/Homo_sapiens.GRCh38.pep.all.fa \
--blastdb $path/VarDB2.0/UniProtApr2018+Ensembl92+GENCODE28.proteins.fasta \
--snpfa $path/VarDB2.0/MSCanProVar_ensemblV79.filtered.fasta \
--genome $path/Homo_sapiens_assembly38.fasta \ #Genome Masked FASTA to BLAT against
--outdir results \ #Output directory
--activationFilter HCD \
--MS2error 0.02 \
--qval 0.01 \
-profile testing --tasks 4 --thread 1\
-resume \ # use it to resume the jobs from the last stopped process.
--PrecursorMassTolerance 10ppm \
--FragmentMethodID 3 \
--mzid2tsvConverter /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MzidToTsvConverter/net45/MzidToTsvConverter.exe
```
The resulting directory will contain the following files 

    H04N32_novel_peptides.gff3
    H04N32_novel_peptides.txt
    H04N32_novel_saav_specai.txt
    H04N32_variant_peptides.txt
    H04N32_variant_precursorError_plot.pdf
    H04N32_variant_specairesult.txt
    ....
    novel_psmtable.txt
    nov_peptidetable.txt
    variant_psmtable.txt
    var_peptidetable.txt

# 5. Score plots 
This Script plots the distribution of retention time,  precuror mass, precursor mass error and match scores of three MS/MS database search tool( MSGF,SpecEValue, Evalue) .
Script :

* `ScorePlotsFromIPAWPipeline.r`

Require files : 

* `novel_psmtable.txt`
# 6.  Add annotation 

Add annotation to novel coding loci (pseudogenes).

Script :

* `IDpick.py`

Require files
* `nov_peptidetable.txt`
* ` mart_export.txt`
```bash
python IDpick.py nov_peptidetable.txt temp.txt sample_nov_peptidetable.IDpick.txt
```
# 7. Count peptides
Count peptides per loci from IDpick file.

Script :

* `CountPeptidesPerLociFromIDpick.r`

Require file:

* `sample_nov_peptidetable.IDpick.txt` 

# 8. Calculate MS1 intensity
we use [moFF](https://github.com/compomics/moFF) published on nature method to extract MS1 intensity(label-free).

1. prepare input files (txt_data) and move them to a subfolder txt_data

   Script :

   * `PrepareTxt_data.r`

   Require file:

   * `novel_psmtable.txt`

2. create a subfolder raw_data and create links to all raw files in this folder.
   If all raw data files are in one folder (without subfolders inside), you can skip this step,using that folder path as the "raw_repo" parameters in the following config file

3. prepare config file
   see example:`config_example.ini`

4. run moFF

   ```bash
   python /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/moFF/moFF/moff_all.py --config_file /path/config_example.in
   ```

# 9. Count reads

The script searchs the reads that support novel peptides from RNA-Seq BAM files.

Script :

- `scam_bams.py`

Require file:

- `novel_peptides.gff3`
- `bam_files_list.txt`  

```bash
python scam_bams.py --input_gff novel_peptides.gff3 --bam_files bam_files_list.txt --output novelpep_readcount.txt
````
