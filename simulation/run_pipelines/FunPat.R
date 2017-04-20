#!/bin/Rscript
## testing R-package FunPat for TC DEG / clustering
############################## library loading ##############################
library(FunPat)
library(edgeR)
###################################################### data input ######################################################
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile1 inFile2 inFileN outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_contr.txt")
    q(save="no")
}

for (idx in 1:(length(args)-1))
{
    outFile <- sub("(contr|treat).txt","FunPat.txt",basename(args[idx]))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- as.numeric(strsplit(par[2],"rep")[[1]])
    TP <- as.numeric(strsplit(par[3],"TP")[[1]])

    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)
    colnames(treat) <- colnames(contr)
    data <- list("contr"=contr,"treat"=treat)
    norm <- lapply(names(data),function(exp){expr=DGEList(data[[exp]]);expr=calcNormFactors(expr);cpm(expr, normalized.lib.sizes=TRUE)})
    names(norm) <- names(data)
    ############################## data pre-processing ##############################
    avg <- lapply(names(norm),function(exp) 
                sapply(1:TP,function(t) 
                    rowMeans(subset(norm[[exp]],select=colnames(norm[[exp]])[((t-1) * repl + 1):((t-1) * repl + repl)]))))
    names(avg) <- names(norm)
    ############################## analysis ##############################
    nC <- avg[["contr"]]
    nD <- avg[["treat"]]
    ## create replicate combinations for treatment and control for each gene and each time point
    replic <- Reduce(rbind,lapply(1:TP,function(idx) cbind(nC[,idx],nD[,idx])))
    
    rank.res <- SEL.TS.AREA(replicates=replic,data1=nC,data2=nD,NAcontrol=TP-1,is.interactive=F,takelog=T)
    outData <- cbind(rank.res$element_ID,rank.res$adjusted_p_value)
    outData <- outData[outData[,2]<=0.01,]
    write.table(outData,file.path(args[length(args)],outFile),row.names=F,col.names=F,quote=F,sep="\t")
}