#!/bin/bash
source /home/zhuyafeng/.bashrc
export NXF_WORK="/ldfssz1/ST_CANCER/CGR/USER/maleyao/IPAW_work/Normal/adrenal_gland/work"
path="/hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/proteogenomics-analysis-workflow/"
path2="/ldfssz1/ST_CANCER/CGR/USER/maleyao/mutationDB"

nextflow run $path/ipaw.hg38.MSGFplus2018.v6.nf --msgfjar /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MSGFplus/MSGFPlus.jar --tdb $path2/VarDB2+Ensembl92.noCOSMIC.revCat.fa --gtf $path/VarDB2.0/VarDB2.gtf --inst 3 --mods labelfree_Mods.txt --mzmldef adrenal_gland_filelist_setnames.txt --knownproteins $path/Ensembl/Homo_sapiens.GRCh38.pep.all.fa --blastdb $path/VarDB2.0/UniProtApr2018+Ensembl92+GENCODE28.proteins.fasta --snpfa $path/VarDB2.0/MSCanProVar_ensemblV79.filtered.fasta --genome $path/Homo_sapiens_assembly38.fasta --outdir results --activationFilter HCD --MS2error 0.02 --qval 0.01 -profile testing --tasks 4 --thread 1 -resume --PrecursorMassTolerance 10ppm --FragmentMethodID 3 --mzid2tsvConverter /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MzidToTsvConverter/net45/MzidToTsvConverter.exe
