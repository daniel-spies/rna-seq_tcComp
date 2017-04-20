#!/bin/Rscript
############ path to directories ############
path <- "~/Documents/phd/data/simulation/study"
outPath <- paste0(path,"/summary_noCutoff")
system(paste0("mkdir -p ",outPath))
methods <- dir(file.path(path,"results_noCutoff"))
files <- sub(paste0(methods[6],".txt"),"",list.files(file.path(path,"results_noCutoff",methods[6]),pattern="*.txt")) # scan files of any of the folders
colHeading <- c("dynb"="score","EBSeqHMM"="score","FunPat"="P","impulseDE"="P","maSigPro"="P","nsgp"="score","pairwise"="P","splineTC"="P","timeSeq"="P")
truth <- read.table(paste0(path,"/stats/SIM_truth.txt"),header=T)
############ combine results for each kind of test ############
for (file in files)
{
    ## list of candidates identified by each method
    results <- list()
    for (type in methods)
    {
        comp <- paste0(path,"/results_noCutoff/",type,"/",file,type,".txt")
        if (!file.exists(comp)) next
        data <- read.table(comp,stringsAsFactors=F,row.names=NULL)
        if (ncol(data) > 2) data[,3:ncol(data)] <- NULL # remove unwanted column as only rownames + P / score in first column needed
        colnames(data) <- c("feature",paste0(type,":",colHeading[[type]]))
        results <- c(results,list(data))
        names(results)[length(results)] <- type
    }
    outData <- Reduce(function(x,y) merge(x,y,by="feature",all=T),results)
    write.table(outData,paste0(outPath,"/",file,"results.txt"),quote=F,row.name=F,sep="\t")
}