#!/bin/Rscript

## set paths and load identifier files
dir <- "~/Documents/phd/data/simulation/study"
files <- list.files(file.path(dir,"results"),pattern=".txt",recursive=T)
files <- files[-grep("(noise|GSE|block|raw)",files)]
DEG <- scan(file.path(dir,"stats/DEG_IDs.txt"),what="",sep="\n")
NEG <- scan(file.path(dir,"stats/NEG_IDs.txt"),what="",sep="\n")
SIM <- read.table(file.path(dir,"stats/SIM_IDs.txt"),header=F,stringsAsFactors=F)

## preparing output matrix
outData <- matrix(ncol=11,nrow=0)
colnames(outData) <- c("TP","FP","FN","TN","sensitivity","specificity","precision","FNR","FDR","accuracy","F1")

for (file in files)
{
    data <- read.table(file.path(dir,"results",file),header=F,quote='""')$V1
    
    ## only TRAP was run on gene names, convert other gene_numbers to gene_symbols
    if(dirname(file) != "TRAP")
    {
        idx <- match(data,SIM$V3)
        data <- SIM[idx,2]
    }
        
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
write.csv(outData,file.path(dir,"stats/simStudy_TC_DEG_statistics.csv"),row.names=TRUE,quote=FALSE)