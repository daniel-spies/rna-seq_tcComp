#!/bin/Rscript
library(ROCR)
path <- "~/Documents/phd/data/simulation/study"
truth <- read.table(file.path(path,"stats/SIM_truth.txt"),header=T,row.names=1)

# loading results table
all <- read.table(file.path(path,"summary/30M_3rep_4TP_results.txt"),header=T,row.names=1)
top5 <- all[,c(1,3,4,6,7)]

data <- list("overlaps"=all, "overlaps_top5"=top5)

for (type in names(data))
{
    newData <- rep(1,nrow(data[[type]]))
    names(newData) <- rownames(data[[type]])
    
    roc <- lapply(1:ncol(data[[type]]), function(overlap){
        # create new results table according to overlap
        overlapIdx <- which(rowSums(data[[type]]==0.01) >= overlap)
        overlapData <- newData
        overlapData[overlapIdx] <- 0.01
    
        pred <- prediction(as.numeric(overlapData == 0.01),truth[names(overlapData),"DEG"])
        perf <- performance(pred,"auc")
        unlist(slot(perf,"y.values"))
    })

    names(roc) <- paste0("overlap_",1:ncol(data[[type]]))
    write.table(as.data.frame(roc),paste0(path,"/stats/AUROC_",type,".txt"),row.names=F)
}