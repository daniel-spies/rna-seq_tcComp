#!/bin/Rscript
library(ggplot2)

path <- "~/Documents/phd/data/simulation/study"
tools <- dir(file.path(path,"results/noise"))
meth <- vector("list",length(tools))
names(meth) <- tools

for (tool in names(meth))
{
    files <- list.files(file.path(path,"results/noise",tool),pattern=paste0(tool,".txt"),full.names=T)
    results <- vector("list",length(files))
    noise <- sub(".*_([0-9.]*)wn_.*","\\1",files)
    names(results) = noise
    
    for (file in files)
    {
        wn <- sub(".*_([0-9.]*)wn_.*","\\1",file)
        data <- read.table(file,quote='"')[,1]
        data <- data[!is.na(data)]
        candidates <- sapply(strsplit(data,"_"),"[[",2)
        total <- length(candidates)
        TP <- sum(as.numeric(candidates)<1201)
        FP <- total - TP
        FN <- 1200 - TP
        TN <- 17303 - FP - FN
        TPR <- TP / ( TP + FN )
        FPR <- FP / (FP + TN)
        results[[wn]] <- c(TPR,FPR)
    }
    outData <- t(data.frame(results))
    colnames(outData) <- c("TPR","FPR")
    outData <- cbind(as.data.frame(outData),tool,noise)
    meth[[tool]] <- outData
}

plotData <- Reduce(rbind,meth)
ggplot(plotData,aes(FPR,TPR,group=tool,shape=noise)) + 
    geom_point(aes(shape=noise,colour=tool),size=5) + 
    geom_vline(data=data.frame(x=c(0.05,0.1)),aes(xintercept=x),colour="red",lty="dotted")
dev.print(cairo_pdf,file.path(path,"plots/results/30M_3rep_4TP_noise.pdf"))