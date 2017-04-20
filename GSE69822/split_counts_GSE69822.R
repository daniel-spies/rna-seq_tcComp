#!/bin/Rscript
library(edgeR)

#extract column names from SraRunTable
dir <- "/cclab_nas/dspies/data/simulated/study/data/GSE69822"
table <- read.csv(file.path(dir,"SraRunTable.txt"),sep="\t")
table <- table[order(table$Run_s),]
sampleNames <- paste(sapply(strsplit(table$genotype_s," "),"[[",1),sapply(strsplit(table$timepoint_s," "),"[[",1),1:3,sep="_")
names(sampleNames) <- table$Run_s

data <- read.table(file.path(dir,"GSE69822.counts"),header=T)
data <- data[,c(1,7:ncol(data))]
rownames(data) <- data$Geneid
data$Geneid <- NULL
colnames(data) <- sampleNames[sub(".*.mapped.(SRR[0-9]*).fastq.*","\\1",colnames(data))]

# filter out non expressed genes
data <- data[rowSums(data)>0,]

# split into PTEN and PIK3CA datasets and filter non expressed genes
wt_idx <- 1:18
pik3ca_idx = treat_idx <- 19:36
pten_idx <- 37:54

PIK3CA <- data[,c(wt_idx,pik3ca_idx)]
PIK3CA <- PIK3CA[rowSums(PIK3CA)>0,]
write.table(PIK3CA[,wt_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pik3ca_contr.txt"),row.names=T,col.names=T,quote=F)
write.table(PIK3CA[,treat_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pik3ca_treat.txt"),row.names=T,col.names=T,quote=F)

PTEN <- data[,c(wt_idx,pten_idx)]
PTEN <- PTEN[rowSums(PTEN)>0,]
write.table(PTEN[,wt_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pten_contr.txt"),row.names=T,col.names=T,quote=F)
write.table(PTEN[,treat_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pten_treat.txt"),row.names=T,col.names=T,quote=F)

# write normalized data for splineTC
group <- rep(1:3,each=18)
d <- DGEList(data,group=group)
d <- calcNormFactors(d)
counts <- cpm(d,normalized.lib.sizes=TRUE)

PIK3CA_norm <- counts[,c(wt_idx,pik3ca_idx)]
PIK3CA_norm <- PIK3CA_norm[rowSums(PIK3CA_norm)>0,]
write.table(PIK3CA_norm[,wt_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pik3ca_norm_contr.txt"),row.names=T,col.names=T,quote=F)
write.table(PIK3CA_norm[,treat_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pik3ca_norm_treat.txt"),row.names=T,col.names=T,quote=F)

PTEN_norm <- counts[,c(wt_idx,pten_idx)]
PTEN_norm <- PTEN_norm[rowSums(PTEN_norm)>0,]
write.table(PTEN_norm[,wt_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pten_norm_contr.txt"),row.names=T,col.names=T,quote=F)
write.table(PTEN_norm[,treat_idx],file.path(dir,"30M_3rep_6TP_GSE69822_pten_norm_treat.txt"),row.names=T,col.names=T,quote=F)

# creating WT - time course by sampling tp 0 for contr matrix
wt_contr <- data[,sample(1:3,18,replace=T)]
wt_contr_norm <- counts[,sample(1:3,18,replace=T)]
colnames(wt_contr) = colnames(wt_contr_norm) <- colnames(data)[wt_idx]

write.table(data[,wt_idx],file.path(dir,"30M_3rep_6TP_GSE69822_wt_treat.txt"),row.names=T,col.names=T,quote=F)
write.table(wt_contr,file.path(dir,"30M_3rep_6TP_GSE69822_wt_contr.txt"),row.names=T,col.names=T,quote=F)

write.table(counts[,wt_idx],file.path(dir,"30M_3rep_6TP_GSE69822_wt_norm_treat.txt"),row.names=T,col.names=T,quote=F)
write.table(wt_contr_norm,file.path(dir,"30M_3rep_6TP_GSE69822_wt_norm_contr.txt"),row.names=T,col.names=T,quote=F)