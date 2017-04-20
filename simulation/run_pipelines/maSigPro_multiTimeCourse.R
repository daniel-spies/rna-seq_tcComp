#!/bin/Rscript
###################################################### library loading ######################################################
library(maSigPro)
###################################################### data input ######################################################
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile1 inFile2 inFileN outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_[contr].txt - supply file name for control!")
    q(save="no")
}

# parameters for estimation and DE/clustering
deg <- 3
Q <- 0.01
meth <- "backward"

for (idx in 1:(length(args)-1))
{
    outFile <- sub("(contr|treat).txt","maSigPro.txt",basename(args[idx]))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- strsplit(par[2],"rep")[[1]]
    TP <- strsplit(par[3],"TP")[[1]]
    timepoints <- paste(rep(((1:TP)-1)*3,each=repl),"h",sep="")
    header <- paste(timepoints,rep(1:repl),sep="-")

    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)
    
    # merge data frames
    data <- merge(contr,treat,by="row.names")
    rownames(data) <- data$Row.names
    data$Row.names <- NULL
    colnames(data) <- paste(header,rep(c("contr","treat"),each=as.numeric(TP)*as.numeric(repl)),sep="_")
   
    # remove non expressed genes
    data <- data[rowSums(data)>0,]

    # create design matrix
    ctrl <- rep(c(1,0),each=length(header))
    mat <- cbind(Time=as.numeric(sub("h","",timepoints)),Replicate=rep(1:(as.integer(TP)*2),each=as.integer(repl)),Control=ctrl,Treatment=as.numeric(ctrl == 0))
    rownames(mat) <- colnames(data)

    # run differential expression analysis
    NBp <- p.vector(data,design=make.design.matrix(mat,degree=deg),counts=TRUE,Q=Q,MT.adjust="BH")
    NBt <- T.fit(NBp,step.method=meth)
    get<-get.siggenes(NBt, vars="groups")
    diff <- get$sig.genes$TreatmentvsControl$sig.pvalues[,c(1,2)]
    write.table(diff,paste(args[length(args)],outFile,sep="/"),row.names=T,col.names=F,quote=F)
}