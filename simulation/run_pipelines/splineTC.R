#!/bin/Rscript
###################################################### library loading ######################################################
library(splineTCDiffExpr)
###################################################### data input ######################################################
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile1 inFile2 inFileN outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_[contr/treat].txt")
    q(save="no")
}

path <- "/Users/dspies/Documents/phd/data/simulation/study"
outPath <- args[length(args)]

for (idx in 1:(length(args)-1))
{
    print(paste("analyzing:",basename(args[idx])))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- as.numeric(strsplit(par[2],"rep")[[1]])
    TP <- as.numeric(strsplit(par[3],"TP")[[1]])
    time <- ((1:TP)-1)*3
    pVal <- 0.01
    header <- paste(rep(((1:TP)-1)*3,each=repl),"h-",rep(1:repl),sep="")

    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)

    inData <- cbind(contr,treat)
    type <- c("contr","treat")
    colnames(inData) <- paste(rep(type,each=length(header)),rep(header,2),sep="_")
    
    design <- data.frame(row.names=colnames(inData),
        "SampleName"=colnames(inData),
        "Time"=rep(rep(time,each=repl),2), 
        "Treatment"=rep(type,each=TP*repl),
        "Replicate"=rep(1:repl,TP*2))
    phenoData <- new("AnnotatedDataFrame",data=design)
    data <- ExpressionSet(assayData=as.matrix(inData),phenoData=phenoData)
    diffExprs <- splineDiffExprs(eSetObject = data, df = TP-1,cutoff.adj.pVal = pVal, reference = type[1],intercept = TRUE)

    outFile <- sub("(contr|treat).txt","splineTC.txt",basename(args[idx]))
    write.table(diffExprs,paste(args[length(args)],outFile,sep="/"),row.names=T,col.names=F,sep="\t")
}