#!/bin/Rscript
## lmms R package for TC DEG 
############################## library loading ##############################
library(lmms)
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
    outFile <- sub("(contr|treat).txt","lmms.txt",basename(args[idx]))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- as.numeric(strsplit(par[2],"rep")[[1]])
    TP <- as.numeric(strsplit(par[3],"TP")[[1]])

    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)
    colnames(treat) <- colnames(contr)

    data <- merge(treat,contr,all=T,by=0)
    rownames(data) <- data$Row.names
    data[,"Row.names"] <- NULL
    
    time <- rep(rep(1:TP,each=repl),2)
    group <- rep(1:2,each=TP*repl)
    sampleID <- sub("y","contr",sub("x","treat",colnames(data)))

    degs <- lmmsDE(t(data),time,sampleID,group,type="time*group",experiment="all",basis="cubic",numCores=25)
                    
    outData <- degs@DE[degs@DE$adj.Group_Time<=0.01,c("Molecule","adj.Group_Time")]
    write.table(outData,file.path(args[length(args)],outFile),row.names=F,col.names=F,quote=F,sep="\t")
}