setwd("D:/project/CPTAC/stomachCancer/S022_pg_nature2014/s22_results/")

df.nov.psm = read.table("novel_psmtable.txt",sep = "\t",header = T,stringsAsFactors = F,
                        quote = "",comment.char = "")

### change pdf file name and the title of each plot                       
pdf("../novel.scoreplot.pdf",width = 9, height = 7, useDingbats = F)
par(mfrow=c(2,3))
hist(df.nov.psm$Retention.time.min.,breaks = 20,main="s22 novel peptides",xlab="PSM retention time (min)")
hist(df.nov.psm$Precursor,breaks = 20,main="s22 novel peptides",xlab="PSM precuror mass")
hist(df.nov.psm$PrecursorError.ppm.,breaks = 20,main="s22 novel peptides", xlab="PSM precursor mass error (ppm)")
hist(df.nov.psm$MSGFScore,breaks = 20,main="s22 novel peptides",xlab="PSM MSGFscore")
hist(-log10(df.nov.psm$SpecEValue),breaks = 20,main="s22 novel peptides",xlab="PSM SpecEValue")
hist(-log10(df.nov.psm$EValue),breaks = 20,main="s22 novel peptides",xlab="PSM Evalue")
dev.off()

df.var.psm = read.table("variant_psmtable.txt",sep = "\t",header = T,stringsAsFactors = F,
                        quote = "",comment.char = "")

### change pdf file name and the title of each plot                      
pdf("../variant.scoreplot.pdf",width = 9, height = 7, useDingbats = F)
par(mfrow=c(2,3))
hist(df.nov.psm$Retention.time.min.,breaks = 20,main="s22 variant peptides",xlab="PSM retention time (min)")
hist(df.nov.psm$Precursor,breaks = 20,main="s22 variant peptides",xlab="PSM precuror mass")
hist(df.nov.psm$PrecursorError.ppm.,breaks = 20,main="s22 variant peptides", xlab="PSM precursor mass error (ppm)")
hist(df.nov.psm$MSGFScore,breaks = 20,main="s22 variant peptides",xlab="PSM MSGFscore")
hist(-log10(df.nov.psm$SpecEValue),breaks = 20,main="s22 variant peptides",xlab="PSM SpecEValue")
hist(-log10(df.nov.psm$EValue),breaks = 20,main="s22 variant peptides",xlab="PSM Evalue")
dev.off()