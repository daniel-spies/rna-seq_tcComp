#!/bin/Rscript
###################################### library loading ######################################
library(DESeq2)
###################################### paths and input ######################################
path <- "/Users/dspies/Documents/phd/data/simulation/study"
pten <- read.table(file.path(path,"run_data/GSE69822/30M_3rep_6TP_GSE69822_pten_treat.txt"),row.names=1)
pik3a <- read.table(file.path(path,"run_data/GSE69822/30M_3rep_6TP_GSE69822_pik3ca_treat.txt"),row.names=1)
wt <- read.table(file.path(path,"run_data/GSE69822/30M_3rep_6TP_GSE69822_wt_treat.txt"),row.names=1)
###################################### parameters ######################################
repl <- 3
TP <- 6
pVal <- 0.01
###################################### analysis ######################################
# comparing WT with PTEN-KO / PIK3A-KO at T=0
wt_pik <- merge(wt[,1:3],pik3a[,1:3],by.x=0,by.y=0,all=T,sort=T)
wt_pik[is.na(wt_pik)] <- 0
rownames(wt_pik) <- wt_pik$Row.names
wt_pik$Row.names <- NULL
wt_pik <- wt_pik[rowSums(wt_pik) >0,]

wt_pten <- merge(wt[,1:3],pten[,1:3],by.x=0,by.y=0,all=T,sort=T)
wt_pten[is.na(wt_pten)] <- 0
rownames(wt_pten) <- wt_pten$Row.names
wt_pten$Row.names <- NULL
wt_pten <- wt_pten[rowSums(wt_pten) >0,]

t0 <- list("PIK3CA"=wt_pik,"PTEN"=wt_pten)

for (type in names(t0))
{
    design <- data.frame(condition=as.factor(rep(c("wt","type"),each=repl)))
    cds <- DESeqDataSetFromMatrix(countData = t0[[type]], colData = design, design = ~condition)
    cds <- DESeq(cds)
    outData <- subset(results(cds),padj <= pVal)
    outFile <- paste0(path,"/results/GSE69822/",type,"_WT_T0.txt")
    write.table(unique(rownames(outData)),outFile,row.names=F,col.names=F,quote=F)
}

# WT TC compared to T=0
# design <- data.frame(time=as.factor(rep(c("0", "15", "40", "90", "180", "300"),each=repl)))
# cds <- DESeqDataSetFromMatrix(countData = wt, colData = design, design = ~time)
# cds <- DESeq(cds)
# cdsLRT <- nbinomLRT(cds, reduced = ~1)
# outData <- subset(results(cdsLRT),padj <= pVal)
# outFile <- paste0(path,"/results/GSE69822/WT_TC/results.txt")
# write.table(unique(rownames(outData)),outFile,row.names=F,col.names=F,quote=F)