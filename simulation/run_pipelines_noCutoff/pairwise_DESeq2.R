#!/usr/lib/R/bin/Rscript
###################################################### library loading ######################################################
library(DESeq2)
library(edgeR)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile1 inFile2 inFileN outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_[contr].txt - supply file name for control!")
    q(save="no")
}
pVal <- 0.01
thresh <- 5
###################################################### edgeR analysis ######################################################

for (idx in 1:(length(args)-1))
{
    outFile <- sub("(contr|treat).txt","pairwise_DESeq2.txt",basename(args[idx]))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- as.numeric(strsplit(par[2],"rep")[[1]])
    TP <- as.numeric(strsplit(par[3],"TP")[[1]])

    timepoints <- paste(rep(1:TP,each=repl),"min",sep="")
    header <- paste(timepoints,rep(1:repl),sep="-")

    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)
    colnames(treat) = colnames(contr) = header

    # remove non expressed genes of each table ( otherwise the defiance function in glm.will throw an error)
    data <- merge(contr,treat,by=0,all=T)
    rownames(data) <- data$Row.names
    data$Row.names <- NULL
    data <- data[rowSums(cpm(data))>thresh,]

    designMatrix <- data.frame(condition=as.factor(rep(c("contr","treat"),each=repl*TP)),time=as.factor(rep(timepoints,2)))
    cds <- DESeqDataSetFromMatrix(countData = data, colData = designMatrix, design = ~condition + time + condition:time)

    # perform the main computation
    cds <- DESeq(cds)
    resultsNames(cds)

    # likelihood ratio test
    cdsLRT <- nbinomLRT(cds, reduced = ~time)
    outData <- subset(results(cdsLRT),padj <= pVal)
    write.table(unique(rownames(outData)),file.path(args[length(args)],outFile),row.names=F,col.names=F,quote=F)
}