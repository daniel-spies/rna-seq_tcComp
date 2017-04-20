############################################# library loading ###########################################
library(edgeR)
library(dplyr)
############################################# data loading ###########################################
cheung <- read.table("~/Documents/phd/data/simulation/study/datasets/cheung_count_table.txt",header=TRUE,stringsAsFactors=FALSE,row.names=1)
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
nSim <- 30000
nRep <- 3
nTP <- 4
multTP <- c(2,3) # multiplicator for time points, higher sampling rate, in theory just more replicates but depicted as time points
################################### simulating genes by drawing from a neg. binomial distribution ####################################
mu <- rnbinom(n=nSim,mu=mean(rowMeans(f$counts)),size=f$common.dispersion)
idx <- sapply(mu,function(x) which.min(abs(pair_f[,1]-min(x,max(pair_f[,1])))))
counts <- data.frame(t(sapply(idx,function(x) replicate(nTP,round(mean(rnbinom(n=min(round(max(pair_f[x,"mu"]/4,1)),500),mu=pair_f[x,"mu"],size=pair_f[x,"d"])))))))
counts <- counts * 4 # increase library depth
data <- counts[sample(1:nrow(counts),20000),] # sampling 20k genes
data <- data[rowSums(data)>0,]
################################### simulating patterns ####################################

# multipliers for patterns
low_up <- c(2,1.2)
low_down <- c(1/2,1/1.2)
high_up <- c(6,3.5)
high_down <- c(1/6,1/3.5)
early_slow_up <- c(2,6,6,6)
early_slow_down <- c(1/2,1/6,1/6,1/6)
mid_slow_up <- c(2,6,6)
mid_slow_down <- c(1/2,1/6,1/6)
mid_fast_up <- c(6,6,6)
mid_fast_down <- c(1/6,1/6,1/6)
grad_lin_up <- c(2,4,6)
grad_lin_down <- c(1/2,1/4,1/6)
grad_cur_up <- c(2,3.5,4.5)
grad_cur_down <- c(1/2,1/3.5,1/4.5)
up_down_slow <- c(6,3,1.2)
up_down_fast <- c(6,1.2)

# row index of DEGs
row <- list(
    up_early_low = seq(1,50),   # selecting transcripts (rows) for DE over 2 timepoints
    up_early_high = seq(51,100),
    up_mid_low = seq(101,150),
    up_mid_high = seq(151,200),
    up_late_low = seq(201,250),
    up_late_high = seq(251,300),
    down_early_low = seq(301,350),
    down_early_high = seq(351,400),
    down_mid_low = seq(401,450),
    down_mid_high = seq(451,500),
    down_late_low = seq(501,550),
    down_late_high = seq(551,600),    
    up_early_slow = seq(601,650),   # selecting transcripts (rows) for DE over 3 timepoints
    up_mid_slow = seq(651,700),
    down_early_slow = seq(701,750),
    down_mid_slow= seq(751,800),
    up_mid_fast = seq(801,850),
    down_mid_fast= seq(851,900),
    up_grad_lin = seq(901,950),
    up_grad_cur = seq(951,1000),
    down_grad_lin = seq(1001,1050),
    down_grad_cur = seq(1051,1100),
    up_down_slow = seq(1101,1150),
    up_down_fast = seq(1151,1200))

# column index of DEGs
col <- list(
    up_early_low = c(1,1),   # selecting transcripts (rows) for DE over 2 timepoints
    up_early_high = c(1,1),
    up_mid_low = c(2,2),
    up_mid_high = c(2,2),
    up_late_low = c(3,3),
    up_late_high = c(3,3),
    down_early_low = c(1,1),
    down_early_high = c(1,1),
    down_mid_low = c(2,2),
    down_mid_high = c(2,2),
    down_late_low = c(3,3),
    down_late_high = c(3,3),
    up_early_slow = c(1,1,1,1),   # selecting transcripts (rows) for DE over 3 timepoints
    up_mid_slow = c(2,2,2),
    down_early_slow = c(1,1,1,1),
    down_mid_slow= c(2,2,2),
    up_mid_fast = c(2,2,2),
    down_mid_fast= c(2,2,2),
    up_grad_lin = c(2,2,2),
    up_grad_cur = c(2,2,2),
    down_grad_lin = c(2,2,2),
    down_grad_cur = c(2,2,2),
    up_down_slow = c(2,2,2),
    up_down_fast = c(2,2))

# pattern modification of each category
mod <- list(
    up_early_low = low_up,   # selecting transcripts (rows) for DE over 2 timepoints
    up_early_high = high_up,
    up_mid_low = low_up,
    up_mid_high = high_up,
    up_late_low = low_up,
    up_late_high = high_up,
    down_early_low = low_down,
    down_early_high = high_down,
    down_mid_low = low_down,
    down_mid_high = high_down,
    down_late_low = low_down,
    down_late_high = high_down,
    up_early_slow = early_slow_up,   # selecting transcripts (rows) for DE over 3 timepoints
    up_mid_slow = mid_slow_up,
    down_early_slow = early_slow_down,
    down_mid_slow= mid_slow_down,
    up_mid_fast = mid_fast_up,
    down_mid_fast= mid_fast_down,
    up_grad_lin = grad_lin_up,
    up_grad_cur = grad_cur_up,
    down_grad_lin = grad_lin_down,
    down_grad_cur = grad_cur_down,
    up_down_slow = up_down_slow,
    up_down_fast = up_down_fast)

## create replicates for more time points
for (TP in 1:length(multTP))
{
    # simulate replicates for control data set
    idx <- sapply(unlist(t(data)),function(x) which.min(abs(pair_f[,"mu"]-x)))
    outData <- matrix(sapply(idx,function(x) replicate(nRep*multTP[TP],round(mean(rnbinom(n=min(round(max(pair_f[x,"mu"]/4,1)),500),mu=pair_f[x,"mu"],size=pair_f[x,"d"]))))),ncol=nRep * nTP * multTP[TP],byrow=TRUE)
    colnames(outData) <- paste("TP",rep(1:(nTP*multTP[TP]),each=nRep),"_",rep(seq(1:nRep),nTP*multTP[TP]),sep="")
    rownames(outData) <- paste("gene_",seq(1,nrow(outData)),sep="")
    write.table(outData,paste("~/Documents/phd/data/simulation/study/countTables/timepoints/30M_3rep_",multTP[TP]*nTP,"TP_contr.txt",sep=""))

    # simulate pattern of treatment data set
    for (sim in names(row))
    {
        mu <- matrix(round(rep(rowMeans(data[row[[sim]],]),each=length(col[[sim]])) * rep(mod[[sim]],length(row[[sim]]))),ncol=length(col[[sim]]),byrow=T)
        cols <- seq(col[[sim]][1],col[[sim]][1]+ncol(mu)-1)
        data[row[[sim]],cols] <- mu
    }
    # simulate replicates for treatment data set
    idx <- sapply(unlist(t(data)),function(x) which.min(abs(pair_f[,"mu"]-x)))
    outData <- matrix(sapply(idx,function(x) replicate(nRep*multTP[TP],round(mean(rnbinom(n=min(round(max(pair_f[x,"mu"]/4,1)),500),mu=pair_f[x,"mu"],size=pair_f[x,"d"]))))),ncol=nRep * nTP * multTP[TP],byrow=TRUE)
    colnames(outData) <- paste("TP",rep(1:(nTP*multTP[TP]),each=nRep),"_",rep(seq(1:nRep),nTP*multTP[TP]),sep="")
    rownames(outData) <- paste("gene_",seq(1,nrow(outData)),sep="")
    write.table(outData,paste("~/Documents/phd/data/simulation/study/countTables/timepoints/30M_3rep_",multTP[TP]*nTP,"TP_treat.txt",sep=""))
}