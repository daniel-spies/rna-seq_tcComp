#!/bin/Rscript
library(ROCR)
library(ggplot2)
path <- "~/Documents/phd/data/simulation/study"
truth <- read.table(file.path(path,"stats/SIM_truth.txt"),header=T,row.names=1)

# loading results table
all <- read.table(file.path(path,"summary/30M_3rep_4TP_results.txt"),header=T,row.names=1)
top5 <- all[,c(4,5,6,8,9)]

data <- list("overlaps"=all, "overlaps_top5"=top5)
dev.new()
for (type in names(data))
{
    newData <- rep(1,nrow(data[[type]]))
    names(newData) <- rownames(data[[type]])
    
    coords <- lapply(1:ncol(data[[type]]), function(overlap){
        # create new results table according to overlap
        overlapIdx <- which(rowSums(data[[type]]==0.01) >= overlap)
        overlapData <- newData
        overlapData[overlapIdx] <- 0.01
    
        pred <- prediction(as.numeric(overlapData == 0.01),truth[names(overlapData),"DEG"])
        coord <- performance(pred,"tpr","fpr")
        c("x"=unlist(slot(coord,"x.values"))[2],"y"=unlist(slot(coord,"y.values"))[2])
    })

    plotCoords <- data.frame(rbind(t(as.data.frame(coords)),matrix(rep(0,ncol(data[[type]])*2),ncol=2),matrix(rep(1,ncol(data[[type]])*2),ncol=2)))
    colnames(plotCoords) <- c("FPR","TPR")
    plotCoords$group <- rep(paste0("overlap_",1:ncol(data[[type]])),3)
    plotCoords <- plotCoords[order(plotCoords[,1],plotCoords[,2]),]
    
    # plotting ROC
    ggplot(plotCoords) + 
        geom_path(aes_string(x = "FPR", y = "TPR", group = "group", color= "group")) + 
        geom_vline(colour="red", linetype="dashed", size = 0.25, xintercept = c(0.01,0.05,0.1))
    dev.print(cairo_pdf,paste0(path,"/plots/results/AUROC_",type,".pdf"),height=10,width=10)
}