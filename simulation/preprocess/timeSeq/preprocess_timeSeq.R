#!/bin/bash
# preprocess simulation study data for timeSeq
# splitting each data set into 10 sets

setwd('/cluster/home04/biol/dspies/timeSeq/data')
files <- list.files(pattern="*contr.txt")
nSplit <- 10

for (file in files)
{
    par <- strsplit(basename(file),"_")[[1]]
    repl <- as.numeric(strsplit(par[2],"rep")[[1]])
    TP <- as.numeric(strsplit(par[3],"TP")[[1]])
    timepoints <- paste(rep(((1:TP)-1)*3,each=repl),"h",sep="")
    header <- paste(timepoints,rep(1:repl),sep="-")

    contr <- read.table(file,header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",file),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)
    colnames(treat) = colnames(contr) = header
    
    ## split data into sets
    start <- as.integer(seq(1,nrow(treat)+1,length.out=nSplit+1))
    outDir <- sub("_contr.txt","",file)
    system(paste("mkdir -p",outDir))
    for(idx in 1:nSplit)
    {
        ##### transform data
        dat_treat <- Reduce(cbind,lapply(1:TP,function(tp) stack(treat[start[idx]:(start[idx+1]-1),((tp-1)*repl + 1):(tp * repl)])$values))
        dat_contr <- Reduce(cbind,lapply(1:TP,function(tp) stack(contr[start[idx]:(start[idx+1]-1),((tp-1)*repl + 1):(tp * repl)])$values))
        dat <- cbind(dat_contr,dat_treat)
        rownames(dat) <- rep(rownames(contr)[start[idx]:(start[idx+1]-1)],repl)
        colnames(dat) <- paste(rep(c("contr","treat"),each=TP),rep(1:TP,2),sep="_")
        dat <- dat[rowSums(dat)>0,]
        outFile <- file.path(outDir,sub("contr",paste0("block_",idx),file))
        print(paste("writing file",outFile))
        write.table(dat,outFile,row.names=T,col.names=T,quote=F)
    }
    print("----------")
}
