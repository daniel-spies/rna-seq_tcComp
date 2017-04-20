#!/bin/Rscript
library(VennDiagram)
library(grid)
library(gridBase)
library(RColorBrewer)

## set paths and load identifier files
dir <- "~/Documents/phd/data/simulation/study"
files <- list.files(file.path(dir,"test"),pattern=".txt")
SIM <- read.table(file.path(dir,"stats/SIM_IDs.txt"),header=F,stringsAsFactors=F)
DEG <- SIM[1:1200,3]
NEG <- SIM[1201:nrow(SIM),3]

## preparing output matrix
outData <- matrix(ncol=11,nrow=0)
colnames(outData) <- c("TP","FP","FN","TN","sensitivity","specificity","precision","FNR","FDR","accuracy","F1")
dataList <- vector("list",length(files))
names(dataList) <- files

for (file in files)
{
    data <- read.table(file.path(dir,"test",file),header=F,quote='""')$V1
    dataList[[file]] <- data
    P <- unique(data)
    TP <- sum(P %in% DEG)
    FP <- length(P) - TP
    FN <- sum(!(DEG %in% P))
    TN <- sum(!(NEG %in% P))

    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    precision <- TP / (TP + FP)
    FNR <- FN / (TP + FN)
    FDR <- FP / (TP + FP)
    accuracy <- (TP + TN) / (TP + FP + TN + FN)
    F1 <- 2 * TP / (2 * TP + FP + FN)
    result <- c(TP,FP,FN,TN,sensitivity,specificity,precision,FNR,FDR,accuracy,F1)
    outData <- rbind(outData,unlist(result))
}
rownames(outData) <- basename(files)
write.csv(outData,file.path(dir,"stats/edgeR_DESeq_pairwise_statistics.csv"),row.names=TRUE,quote=FALSE)

### plotting Venn diagram
names(dataList) <- 1:length(dataList)
category <- sub('30M_3rep_4TP_(.*).txt','\\1',files)
col <- c(brewer.pal(max(3,length(category)),"YlOrRd"))[1:length(category)]
colStr <- paste("fill=c(\"",paste(paste(col),collapse="\",\""),"\")",sep="")

combs <- unlist(lapply(1:length(dataList), function(j) combn(names(dataList), j, simplify = FALSE)),recursive = FALSE)
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
elements <- lapply(combs, function(i) Reduce(intersect,dataList[i]))
n.elements <- sapply(elements, length)

## graphical parameters
labels <- paste(category,paste("(n=",n.elements[1:length(category)],")",sep=""),sep="\n")
catStr <- paste("category=c(\"",paste(labels,collapse="\",\""),"\")",sep="")
areas <- as.list(unname(n.elements))
areaStr <- paste(areas,collapse=",")
graphStr <- paste("ind=FALSE,cat.cex=rep(2,times=",length(category),"),cex=rep(1.5,times=",length(elements),"),margin=0.1",sep="")

## join parameters into one string for evaluation and plotting
param <- paste(areaStr,catStr,colStr,graphStr,sep=",")
draw <- switch(length(category),"single","pairwise","triple","quad","quintuple")
grid.draw(eval(parse(text=paste0("draw.",draw,".venn(",param,")"))))
grid.text("Comparison of pairwise DEG tests",y=unit(0.05,"npc"),gp=gpar(fontsize=30))
dev.print(device=cairo_pdf,file.path(dir,"stats/edgeR_DESeq_pairwise_venn.pdf"),width=20,height=20)