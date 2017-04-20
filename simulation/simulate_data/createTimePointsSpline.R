############################################# library loading ###########################################
library(lmms)
library(edgeR)
library(DESeq2)
############################################# functions ###########################################
## check row means of contr / treat 
## check values of spline and plot if it makes sense or introduces artifacts / errors

sampleReplicates <- function(mu, nRep, numTP){
    idx <- sapply(mu,function(x) which.min(abs(pair_f[,"mu"]-x)))
    outData <- matrix(sapply(idx,function(x) replicate(nRep,round(mean(rnbinom(n=min(round(max(pair_f[x,"mu"]/4,1)),500),mu=pair_f[x,"mu"],size=pair_f[x,"d"]))))),ncol=nRep * numTP,byrow=TRUE)
    colnames(outData) <- paste("TP",rep(1:numTP,each=nRep),"_",rep(seq(1:nRep),numTP),sep="")
    rownames(outData) <- paste("gene_",seq(1,nrow(outData)),sep="")
    return(outData)
}
############################################# data loading ###########################################
setwd("/cclab_nas/dspies/data/simulated/study")
counts <- list(treat=read.table("run_data/30M_3rep_4TP_treat.txt"),contr=read.table("run_data/30M_3rep_4TP_contr.txt"))
cheung <- read.table("cheung_count_table.txt",header=TRUE,stringsAsFactors=FALSE,row.names=1)
cheung <- cheung[rowSums(cheung)>0,]
############################################# estimating dispersions ###########################################
f <- DGEList(cheung)
f <- estimateGLMCommonDisp(f) 
f <- estimateGLMTrendedDisp(f)
f <- estimateGLMTagwiseDisp(f)
############################################# create mu / disp dictionary ###########################################
pair_f <- data.frame(mu=round(rowMeans(f$counts)),d=f$tagwise.dispersion)
pair_f <- pair_f[order(pair_f[,1],decreasing=TRUE),]
pair_f <- aggregate(.~mu,data=pair_f,FUN=mean)
############################################# simulation parameters ###########################################
nDEG <- 1200
nRep <- 3
nTP <- 4
multTP <- c(2,3)
time <- rep(1:nTP,each=nRep)
time.regular <- lapply(multTP, function(nOut) seq(min(time), max(time), length.out = nTP * nOut))
means <- round(rowMeans(counts[["contr"]][,1:3]))
################################### simulating spline and interpolate additional time points ####################################
dev.new(width=10,height=10)
par(mfrow=c(2,2))

for (type in names(counts))
{
    outData <- switch(type,
        "contr" = { # sample flat line replicates for control samples from means of the first TP only
            print("sampling from TP1 for control samples")
            lapply(time.regular, function(timeR){
                print(paste0(".... sampling ",length(timeR)," TPs"))
                mu <- rep(means, each=length(timeR))
                repData <- sampleReplicates(mu, nRep, length(timeR))
                write.table(repData,paste0("run_data/30M_3rep_",length(timeR),"TP_",type,"_spline.txt"))
                return(repData)
            })
        },
        "treat" = { # fit spline and sample for spline for first nDEG genes and use flat profile for non DEGs
            print("sampling from spline")
            spline <- lmmSpline(data = t(counts[[type]][1:nDEG,]), time = time, sampleID = 1:ncol(counts[[type]]), keepModels = TRUE, numCores = 10)
            nonDEG <- means[(nDEG+1):nrow(counts[[type]])]
            lapply(time.regular, function(timeR){
                print(paste0(".... sampling ",length(timeR)," TPs"))
                data.interpolate <- t(predict(spline, timePredict = timeR, numCores = 10)) # interpolate to desired numbe of time points
                mu <- matrix(rep(nonDEG, each=length(timeR)),nrow=length(timeR),byrow=F)
                data <- cbind(data.interpolate,mu)
                splineData <- sampleReplicates(data, nRep, length(timeR)) # sample replicates
                write.table(splineData,paste0("run_data/30M_3rep_",length(timeR),"TP_",type,"_spline.txt"))
                return(splineData)
            })
        }
    )
    ################################### check dispersion ####################################
    plotList <- lapply(outData, function(counts){
        TP <- ncol(counts) / nRep
        design <- data.frame(tp=as.factor(rep(1:TP,each=nRep)),repl=as.factor(rep(1:nRep,TP)))
        dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~ tp + repl)
        dds <- estimateSizeFactors(dds)
        estimateDispersions(dds)    
    })
    names(plotList) <- paste(type,c(8,12),"TPs",sep="_")
    for (exp in names(plotList))
    {
        plotDispEsts(plotList[[exp]])
        title(main=exp)
    }
}
dev.print(dev=cairo_pdf,"dispersion_spline_samples.pdf")