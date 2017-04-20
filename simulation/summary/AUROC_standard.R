#!/bin/Rscript
library(iCOBRA)
library(ROCS)
library(ggplot2)
path <- "~/Documents/phd/data/simulation/study"

thresh <- c(0.01,0.05,0.1)

# load data into iCOBRA 
t <- COBRAData_from_text(file.path(path,"stats/SIM_truth.txt"), file.path(path,"summary_noCutoff/30M_3rep_4TP_results.txt"), "feature")
t <- calculate_adjp(t, method = "BH")
cobraperf <- calculate_performance(t, binary_truth = "DEG")
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", facetted = TRUE)

# plot ROC
ob <- plot_roc(cobraplot,linewidth=0.5)

plotData <- cbind(cobraperf@fdrtprcurve, FPR=cobraperf@roc$FPR)
points <- as.data.frame(Reduce(rbind,lapply(unique(cobraperf@fdrtprcurve$method), function(tool)
    t(sapply(thresh, function(thr) {
        idx <- min(which(plotData$method == tool & plotData$FDR >= thr))
        plotData[idx, c("TPR","FPR","method")]
    })))))

# format columns for plotting
points$TPR <- as.numeric(points$TPR)
points$FPR <- as.numeric(points$FPR)
points$method <- as.character(points$method)

# add FDR points to ROC plot
ob <- plot_roc(cobraplot,linewidth=0.5)
ob + geom_point(data=points,mapping=aes(x=FPR,y=TPR,colour=method,group=method),size=2)
dev.print(cairo_pdf,file.path(path,"plots/results/30M_3rep_4TP_roc-data_new.pdf"))

# plot TPR / FDR
plot_fdrtprcurve(cobraplot,linewidth=0.5,pointsize=1)
dev.print(cairo_pdf,file.path(path,"plots/results/30M_3rep_4TP_tpr_fdr-data_new.pdf"))

# calculate AUC for full curve
auc <- sapply(c(thresh,seq(0.2,1,0.1)), function(thr)
    sapply(unique(cobraperf@fdrtprcurve$method), function(tool){
        plotData <- cbind(subset(cobraperf@fdrtprcurve, cobraperf@fdrtprcurve$method == tool, c("FDR","TPR","TP","FP")),
                            subset(cobraperf@roc, cobraperf@roc$method == tool, c("FPR")))
        fcauc.fptp(FP=plotData$FPR, TP=plotData$TPR, TDR=1 - plotData$FDR, FDR.cut = thr, do.plot=F)
    }))
colnames(auc) <- paste0("FDR_",c(thresh,seq(0.2,1,0.1)))
write.table(t(round(auc,2)),file.path(path,"stats/AUROC_standard.txt"),row.names=T,col.names=T,quote=F)