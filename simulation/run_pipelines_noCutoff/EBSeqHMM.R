#!/bin/Rscript
###################################################### library loading ######################################################
library(EBSeqHMM)
###################################################### data input ######################################################
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile1 inFile2 inFileN outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_contr.txt")
    q(save="no")
}

# parameters for estimation and DE/clustering
rounds <- 10
FDR <- 1
qTrm <- .25
qCut <- 0
co = 0

for (idx in 1:(length(args)-1))
{
    outFile <- sub("(contr|treat).txt","EBSeqHMM.txt",basename(args[idx]))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- strsplit(par[2],"rep")[[1]]
    TP <- strsplit(par[3],"TP")[[1]]
    time <- factor(rep(1:TP,each=repl))
    header <- paste(rep(((1:TP)-1)*3,each=repl),"h-",rep(1:repl),sep="")

    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)
    colnames(treat) = colnames(contr) = header

    type <- list("contr"=contr,"treat"=treat)
    genes = degs <- list("contr"=c(),"treat"=c())

    for (meth in names(type))
    {
        colnames(type[[meth]]) <- header
        Sizes <- MedianNorm(type[[meth]])
        GeneNormData <- GetNormalizedMat(type[[meth]],Sizes)
        EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneNormData, sizeFactors=Sizes, Conditions=time, UpdateRd=rounds, Qtrm=qTrm, QtrmCut=qCut)
        GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut, FDR=FDR, cutoff=co, OnlyDynamic=FALSE)
        genes[[meth]] <- GeneConfCalls$Overall
    }
    
    all <- as.data.frame(Reduce(rbind,genes))
    all$gene <- rownames(all)
    outData <- aggregate( Max_PP ~ gene, all, max)
    outData <- outData[rev(order(outData$Max_PP)),]
    write.table(outData[,2:3],paste(args[length(args)],outFile,sep="/"),row.names=F,col.names=F,sep="\t",quote=F)
}