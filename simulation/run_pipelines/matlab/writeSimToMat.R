#!/bin/Rscript
############################## libraries ##############################
library(R.matlab)
library(DESeq2)
############################## read input ##############################
#dir=~/Documents/phd/data/simulation/study/countTables
#out=~/Documents/phd/data/simulation/study/matlabTables

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) args <- c("--help")
if("--help" %in% args) cat("Arguments:\n--args inDir outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_[contr/treat].tx")
    
dir <- args[1]
outDir <- args[2]
for (folder in grep("^(?!raw)",list.files(dir),perl=TRUE,value=TRUE))
{
    workDir <- paste(dir,folder,sep="/")
    print(paste("formating read counts:",workDir))
    ctrl <- list.files(workDir,pattern="contr")
    treat <- list.files(workDir,pattern="treat")
    for (idx in seq(1,length(ctrl)))
    {
        ctrl_count <- read.table(paste(dir,folder,ctrl[idx],sep="/"))
        treat_count <- read.table(paste(dir,folder,treat[idx],sep="/"))
        repl <- as.integer(gsub(".*_(.*)rep_.*",'\\1',ctrl[idx]))
        TP <- as.integer(gsub('.*_(.*)TP_.*.txt','\\1',ctrl[[idx]]))
    
        # factors / conditions for size factor estimation
        subject <- colnames(ctrl_count)
        fact <- factor(sub("_[1-9]*","",subject))
        cond <- data.frame(row.names=subject,conditions=fact)

        # set dimensions for matlab array
        len <- nrow(ctrl_count)
        expMat <- vector("list", len)
        names(expMat) <- rownames(ctrl_count)
        conMat <- expMat

        for (index in 1:len)
        {
            expMat[[index]] <- sapply(1:TP,function(col)treat_count[index,((col-1)*repl+1):(repl*col)])
            conMat[[index]] <- sapply(1:TP,function(col)ctrl_count[index,((col-1)*repl+1):(repl*col)])
        }

        # size factor estimation
        tmp_ctrl <- matrix(as.integer(unlist(ctrl_count)),ncol=ncol(ctrl_count),byrow=TRUE)
        tmp_counts <- matrix(as.integer(unlist(treat_count)),ncol=ncol(treat_count),byrow=TRUE)
        dds_ctrl <- DESeqDataSetFromMatrix(countData = tmp_ctrl,colData = cond,design = ~ conditions)
        sf_ctrl <- estimateSizeFactors(dds_ctrl)
        conSize <- matrix(sf_ctrl$sizeFactor,ncol=TP,byrow=TRUE)
        dds_counts <- DESeqDataSetFromMatrix(countData = tmp_counts,colData = cond,design = ~ conditions)
        sf_counts <- estimateSizeFactors(dds_counts)
        treSize <- matrix(sf_counts$sizeFactor,ncol=TP,byrow=TRUE)
        outFile <- gsub("_contr.txt",".mat",ctrl[idx])
        writeMat(paste(outDir,folder,outFile,sep="/"),control=array(unlist(conMat),dim=c(repl,TP,len)),treatment=array(unlist(expMat),dim=c(repl,TP,len)),treSize=treSize,conSize=conSize,labels=rownames(ctrl_count))
    }
}