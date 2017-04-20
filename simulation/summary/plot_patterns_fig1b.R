#!/bin/Rscript
library(ggplot2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,byRow=F) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols),byrow=byRow)
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
path <- "/Users/dspies/Documents/phd/data/simulation/study"
TP <- 4
repl <- 3
colIdx <- lapply(1:TP,function(x) (((x-1)*repl+1):(x*repl)))
rawData <- read.table(file.path(path,"run_data/30M_3rep_4TP_treat.txt"))

patterns <- list(
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

plots <- vector("list",length(patterns))
names(plots) <- names(patterns)
for (cat in names(patterns))
{
    rowIdx <- patterns[[cat]]
    m <- Reduce(rbind,sapply(colIdx,function(col) rowMeans(rawData[rowIdx,col])))
    time <- rep(c(0,3,6,9),each=length(rowIdx))
    group <- rep(1:50,4)
    plotData <- as.data.frame(cbind(m,time,group))
    colnames(plotData) <- c("mean","time","group")
    
    plots[[cat]] <- ggplot(plotData,aes(y=mean,x=time,color=group,group=group)) + 
        geom_line() + 
        theme_bw() +
        theme(panel.background = element_blank(),panel.border = element_blank(),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             legend.position = "none",
             axis.line = element_blank(),
             axis.text.y = element_blank(), axis.text.x = element_blank(),
             axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
             axis.title.x = element_blank(), axis.title.y = element_blank())
}

multiplot(plotlist=plots,cols=4,byRow=TRUE)
dev.print(device=cairo_pdf,paste(path,"plots/preStart/DEG_standard_patterns.pdf",sep="/"),height=40,width=30)