setwd("E:\\protein_3.0\\normal\\IDpick\\")  
dir = list.files(pattern = "*.IDpick.txt")
n = length(dir)
sum_novel_coding_loci<-data.frame()
df_sum_LORF1_2_HUMAN<-data.frame()
df_sum_repeat<-data.frame()
df_mergeIDpick<-data.frame()
df_source<-data.frame()
sum_source<-data.frame()
library("plyr")
for (i in 1:n){
  dfidpick<-read.table(dir[i],sep="\t",header=T,fill=T,stringsAsFactors = F,quote = "",comment.char = "" )
  
  ####qval contronl  qval<=0.01
  df.qval = dfidpick[,grep(".*q.value",colnames(dfidpick),perl = T)]
  df.qval = df.qval[,grep("^(?!.*PSM).*q.value",colnames(df.qval),perl = T)]
  df.qval[df.qval>0.01] = NA
  idex = rowSums(is.na(df.qval))< length(colnames(df.qval))
  #idex = !is.na(df.qval)
  dfidpick = dfidpick[idex,]
  dfidpick[,colnames(df.qval)][dfidpick[,colnames(df.qval)]>0.01 ]= NA
  index <- duplicated(dfidpick[,c("Sequence","Protein")])
  dfuq<-dfidpick[!index,]
  
  ####noncoding
  #df_noncoding<-dfuq[grep(pattern = 'UTR5|UTR3|exonic',invert = TRUE,dfuq$anovar_category),]
  #df_nc_exonic<-dfuq[grep(pattern = 'ncRNA_exonic',dfuq$anovar_category),]
  #dfuq<-rbind.data.frame(df_noncoding,df_nc_exonic)
  
  ###filteridpick Peptides Frequence >=2
  dfFreq_2<-dfuq[which(dfuq$Protein %in% names(which(table(dfuq$Protein)>=2))),]
  write.table(dfFreq_2,file=paste('E:\\protein_3.0\\normal\\filterIDpick2\\filter2.0_',dir[i],sep=''),sep="\t",row.names = F,quote=F)
  
  ###merge
  df_name=data.frame(c(rep(dir[i],nrow(dfFreq_2))))
  dfsingle_Freq_2<-cbind(df_name,dfFreq_2)
  df_mergeIDpick<-rbind.fill(df_mergeIDpick,dfsingle_Freq_2)
  
  ###novel coding
  t = as.data.frame(table(dfuq$Protein))
  n1 = nrow(t[t$Freq==1,])
  n2 = nrow(t[t$Freq>=2 & t$Freq<=3,])
  n3 = nrow(t[t$Freq>3,])
  
  ###function
  n4 = sum(grepl("actin",fixed = T,dfuq$Description))
  n5 = sum(grepl("eukaryotic translation",fixed = T,dfuq$Description))
  n6 = sum(grepl("glyceraldehyde",fixed = T,dfuq$Description))
  n7 = sum(grepl("heat shock",fixed = T,dfuq$Description))
  n8 = sum(grepl("heterogeneous nuclear ribonucleoprotein",fixed = T,dfuq$Description))
  n9 = sum(grepl("immunoglobulin",fixed = T,dfuq$Description))
  n10 = sum(grepl("keratin",fixed = T,dfuq$Description))
  n11= sum(grepl("peptidylprolyl isomerase A",fixed = T,dfuq$Description))
  n12= sum(grepl("ribosomal protein",fixed = T,dfuq$Description))
  n13 = sum(grepl("tubulin",fixed = T,dfuq$Description))
  
  ###LORF1_2_HUMAN
  df_LORF1_2_HUMAN<-dfFreq_2[grep(pattern = 'LORF1_HUMAN|LORF2_HUMAN',dfFreq_2$Description),]
  n_LORF1_HUMAN = sum(grepl("LORF1_HUMAN",fixed = T,dfFreq_2$Description))
  n_LORF2_HUMAN = sum(grepl("LORF2_HUMAN",fixed = T,dfFreq_2$Description))
  df_name=data.frame(c(rep(dir[i],n_LORF1_HUMAN+n_LORF2_HUMAN)))
  df_single_LORF1_2_HUMAN<-cbind(df_name,df_LORF1_2_HUMAN)
  df_sum_LORF1_2_HUMAN<-rbind.fill(df_sum_LORF1_2_HUMAN,df_single_LORF1_2_HUMAN)
  
  ###repeat
  df_repeat<-dfFreq_2[grep(pattern = 'multiple',dfFreq_2$blat_category),]
  n_repeat=sum(grepl("multiple",fixed = T,dfFreq_2$blat_category))
  df_name=data.frame(c(rep(dir[i],n_repeat)))
  df_single_repeat<-cbind(df_name,df_repeat)
  df_sum_repeat<-rbind.fill(df_sum_repeat,df_single_repeat)
  ###
  novel_coding_loci<-data.frame(datasets=dir[i],one_peptide = n1,two_peptides=n2,three_peptide3=n3,
                                actin=n4,eukaryotic_translation_elongation_factor=n5,	GAPDH=n6,
                                heat_shock_protein=n7,heterogeneous_nuclear_ribonucleoprotein=n8,immunoglobulin=n9,
                                keratin=n10,peptidylprolyl_isomerase_A=n11,ribosomal_protein=n12,tublin=n13,
                                LORF1_HUMAN=n_LORF1_HUMAN,LORF2_HUMAN=n_LORF2_HUMAN,repea=n_repeat)
  sum_novel_coding_loci<-rbind(sum_novel_coding_loci,novel_coding_loci)
  
  ###source
  n21=nrow(dfFreq_2[dfFreq_2$anovar_category=="downstream",])				
  n22 =nrow(dfFreq_2[dfFreq_2$anovar_category=="exonic",])
  n23 =nrow(dfFreq_2[dfFreq_2$anovar_category=="exonic-intronic",])
  n24 =nrow(dfFreq_2[dfFreq_2$anovar_category=="exonic-ncRNA_exonic",])
  n25 =nrow(dfFreq_2[dfFreq_2$anovar_category=="exonic-upstream",])
  n26 =nrow(dfFreq_2[dfFreq_2$anovar_category=="exonic-UTR5",])
  n27 =nrow(dfFreq_2[dfFreq_2$anovar_category=="intergenic",])
  n28 =nrow(dfFreq_2[dfFreq_2$anovar_category=="intergenic-upstream",])
  n29 =nrow(dfFreq_2[dfFreq_2$anovar_category=="intronic",])
  n30=nrow(dfFreq_2[dfFreq_2$anovar_category=="intronic-exonic",])
  n31=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_exonic",])
  n32=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_exonic-intergenic",])
  n33=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_exonic-ncRNA_intronic",])
  n34=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_exonic-upstream",])
  n35=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_intronic",])
  n36=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_intronic-ncRNA_exonic",])
  n37=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_splicing-ncRNA_intronic",])
  n38=nrow(dfFreq_2[dfFreq_2$anovar_category=="upstream",])
  n39=nrow(dfFreq_2[dfFreq_2$anovar_category=="upstream-intergenic",])
  n40=nrow(dfFreq_2[dfFreq_2$anovar_category=="UTR3",])
  n41=nrow(dfFreq_2[dfFreq_2$anovar_category=="UTR5",])
  n42=n41=nrow(dfFreq_2[dfFreq_2$anovar_category=="ncRNA_exonic-exonic",])
  n43=n41=nrow(dfFreq_2[dfFreq_2$anovar_category=="UTR5-intronic",])
  pseudogene=sum(grepl("pseudogene",fixed = T,dfFreq_2$Description))
  df_source<-data.frame(datasets=dir[i],downstream=n21,exonic=n22,'exonic-intronic'=n23,'exonic-ncRNA_exonic'=n24,'exonic-upstream'=n25,'exonic-UTR5'=n26,intergenic=n27,'intergenic-upstream'=n28,intronic=n29,'intronic-exonic'=n30,
                        'ncRNA_exonic'=n31,'ncRNA_exonic-intergenic'=n32,'ncRNA_exonic-ncRNA_intronic'=n33,
                        'ncRNA_exonic-upstream'=n34,ncRNA_intronic=n35,'ncRNA_intronic-ncRNA_exonic'=n36,'ncRNA_splicing-ncRNA_intronic'=n37,upstream=n38,'upstream-intergenic'=n39,UTR3=n40,UTR5=n41,
                        'ncRNA_exonic-exonic'=n42,'UTR5-intronic'=n43,pseudogene=pseudogene,LORF1_HUMAN=n_LORF1_HUMAN,LORF2_HUMAN=n_LORF2_HUMAN)
  sum_source<-rbind(sum_source,df_source)
  
}

write.table(sum_novel_coding_loci,"sum_novel_fun_coding_loci.txt",sep="\t",row.names = F,quote = F)
write.table(df_sum_LORF1_2_HUMAN,"sum_LORF1_2_HUMAN.txt",sep="\t",row.names = F,quote = F)
write.table(df_sum_repeat,"sum_repeat.txt",sep="\t",row.names = F,quote = F)
write.table(df_mergeIDpick,"mergeIDpick.txt",sep="\t",row.names = F,quote = F)
write.table(sum_source,"source.txt",sep="\t",row.names = F,quote = F)

new<-df_mergeIDpick[,1:44]
write.table(new,"mergeIDpick.txt",sep="\t",row.names = F,quote = F)
