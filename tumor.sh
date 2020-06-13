#!/bin/bash

export NXF_WORK="/hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/IPAW_work_folder/s37_PNNL/work_s37_PNNL"
path="/hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/proteogenomics-analysis-workflow/"

nextflow run $path/ipaw.hg38.MSGFplus2018.nf --msgfjar /hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MSGFplus/MSGFPlus.jar --tdb $path/VarDB2.0/VarDB2+coadreadMut+Ensembl92.revCat.fa --gtf $path/VarDB2.0/VarDB2.gtf --inst 3 --mods tmt_Mods.txt --mzmldef s37_PNNL_filelist_setnames.txt --knownproteins $path/Ensembl/Homo_sapiens.GRCh38.pep.all.fa --blastdb $path/VarDB2.0/UniProtApr2018+Ensembl92+GENCODE28.proteins.fasta --snpfa $path/VarDB2.0/MSCanProVar_ensemblV79.filtered.fasta --genome $path/Homo_sapiens_assembly38.fasta --isobaric tmt10plex --denoms 'set01:131 set02:131 set03:131 set04:131 set05:131 set06:131 set07:131 set08:131 set09:131 set10:131 set11:131 set12:131 set13:131 set14:131 set15:131 set16:131 set17:131 set18:131 set19:131 set20:131 set21:131 set22:131' --outdir s37_PNNL_results -profile testing
