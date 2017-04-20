#!/usr/lib/R/bin/Rscript
###################################################### library loading ######################################################
library(edgeR)
###################################################### data input  and parameters######################################################
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile1 inFile2 inFileN outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_[contr].txt - supply file name for control!")
    q(save="no")
}
###################################################### edgeR analysis ######################################################

for (idx in 1:(length(args)-1))
{
    outFile <- sub("(contr|treat).txt","pairwise.txt",basename(args[idx]))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- as.numeric(strsplit(par[2],"rep")[[1]])
    TP <- as.numeric(strsplit(par[3],"TP")[[1]])

    timepoints <- paste(rep(((1:TP)-1)*3,each=repl),"h",sep="")
    header <- paste(timepoints,rep(1:repl),sep="-")

    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)
    colnames(treat) = colnames(contr) = header

    # remove non expressed genes of each table ( otherwise the defiance function in glm.will throw an error)
    data <- merge(contr,treat,by=0,all=T)
    rownames(data) <- data$Row.names
    data$Row.names <- NULL
    DEG <- vector("list",TP)
    names(DEG) <- paste("TP",1:TP,sep="_")

    for (t in 0:(TP-1))
    {
        contr_idx <- c(((t*repl) + 1):((t+1)*repl))
        treat_idx <- c(((TP*repl) + repl*t + 1):((TP*repl) + repl*(t+1)))
        counts <- data[,c(contr_idx,treat_idx)]
        design<- data.frame(condition=as.factor(rep(c("contr","treat"),each=repl)))
        model <- model.matrix(~condition,data=design)

        d <- DGEList(counts,group=c(rep(1,repl),rep(2,repl)))
        d <- calcNormFactors(d)
        d <- estimateGLMCommonDisp(d)
        d <- estimateGLMTrendedDisp(d)
        d <- estimateGLMTagwiseDisp(d)
        et <- exactTest(d)
        table <- et$table[order(et$table$PValue),]
        table$gene <- rownames(table)
        DEG[[t+1]] <- table
    }
    outData <- aggregate(PValue ~ gene, Reduce(rbind,DEG), min)
    outData <- outData[order(outData$PValue),]
    write.table(outData,file.path(args[length(args)],outFile),row.names=F,col.names=F,quote=F)
}