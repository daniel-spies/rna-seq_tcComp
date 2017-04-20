#!/bin/Rscript
######################################################
library(timeSeq)
###################################################### data input ######################################################
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile \n\nno file supplied!")
    q(save="no")
}

print("setting parametes")
pVal <- 0.01
TP <- as.integer(sub(".*_([0-9]*)TP_.*","\\1",args[1]))
outFile <- sub(".txt","_timeSeq.txt",basename(args[1]))
outDir <- sub("data","results",dirname(args[1]))
system(paste("mkdir -p",outDir))

print(paste("reading data",args[1]))
input <- read.table(args[1],header=T,stringsAsFactors=F,row.names=NULL)
data <- as.matrix(input[,2:ncol(input)],ncol=TP*2)
rownames(data) <- input$row.names

print("DEG analysis")
model.fit <- timeSeq(data,rep(1:2,each=TP),rownames(data),reads=NULL,exon.level=F,gene.level=T,p.values=T,n_cores=48)
outData <-  model.fit$pvalue[model.fit$pvalue$NPDE <= pVal | model.fit$pvalue$PDE <= pVal,]

print("writing results")
write.table(outData,file.path(outDir,outFile),row.names=T,col.names=T,quote=F)