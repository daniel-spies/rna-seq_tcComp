#!/bin/Rscript

############ path to directories ############
path <- "~/Documents/phd/data/simulation/study"
outPath <- paste0(path,"/summary")
system(paste0("mkdir -p ",outPath))
methods <- dir(file.path(path,"results"))[-c(4,6)] # exclude noise and GSE data set
files <- sub(paste0(methods[3],".txt"),"",list.files(file.path(path,"results",methods[3])))
############ writing ground truth file ############
## include category for fractioning
ids <- read.table(paste0(path,"/stats/SIM_IDs.txt"),stringsAsFactors=F)[3]
table <- cbind(ids,c(rep(1,1200),rep(0,nrow(ids)-1200)))

cat_time <- c(
    rep(c(rep("early",100),rep("mid",100),rep("late",100)),2), rep(c(rep("early",50),rep("mid",50)),2), rep("mid",100), rep("grad",200), rep("mixed",100),rep("non-DEG",nrow(ids)-1200))
cat_type <- c(rep(c(rep("low",50),rep("high",50)),6), rep("slow",200),rep("fast",100), rep("grad",200), rep("slow",50),rep("fast",50),rep("non-DEG",nrow(ids)-1200))
truth <- cbind(table,cat_time,cat_type)
colnames(truth) <- c("feature","DEG","cat_time","cat_type")
write.table(truth,paste0(path,"/stats/SIM_truth.txt"),quote=F,row.name=F,sep="\t")
############ combine results for each kind of test ############
for (file in files)
{
    ## list of candidates identified by each method
    results <- list()
    for (type in methods)
    {
        comp <- paste0(path,"/results/",type,"/",file,type,".txt")
        if (!file.exists(comp)) next
        data <- cbind(read.table(comp,stringsAsFactors=F)[1],0.01)
        colnames(data) <- c("feature",paste0(type,":adjP"))
        results <- c(results,list(data))
        names(results)[length(results)] <- type
    }
    outData <- Reduce(function(x,y) merge(x,y,by="feature",all=T),results)
    ## add Negativ IDs for complete graphs
    N_IDs <- setdiff(truth$feature,outData$feature)
    N <- cbind(N_IDs,data.frame((matrix(NA,nrow=length(N_IDs),ncol=ncol(outData)-1))))
    colnames(N) <- colnames(outData)
    outData <- rbind(outData,N)
    outData[is.na(outData)] <- "1.00"
    write.table(outData,paste0(outPath,"/",file,"results.txt"),quote=F,row.name=F,sep="\t")
}