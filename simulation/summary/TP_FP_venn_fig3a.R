#!/bin/Rscript
library(VennDiagram)
library(grid)
library(gridBase)
library(RColorBrewer)

## set paths and load identifier files
dir <- "~/Documents/phd/data/simulation/study"
files <- list.files(file.path(dir,"results"),pattern="30M_3rep_4TP.*txt",recursive=T)
files <- files[-grep("(noise|GSE|EBSeqHMM|nsgp|timeSeq)",files)]
DEG <- scan(file.path(dir,"stats/DEG_IDs.txt"),what="",sep="\n")
NEG <- scan(file.path(dir,"stats/NEG_IDs.txt"),what="",sep="\n")
SIM <- read.table(file.path(dir,"stats/SIM_IDs.txt"),header=F,stringsAsFactors=F)

data <- lapply(files, function(file)
{
    data <- read.table(file.path(dir,"results",file),header=F,quote='""')$V1
    idx <- match(data,SIM$V3)
    data <- SIM[idx,2]
    
    P <- unique(data)
    TP <- P[P %in% DEG]
    FP <- P[P %in% NEG]
    result <- list("TP"=TP,"FP"=FP)
    return(result)
})
names(data) <- sub('.*/30M_3rep_4TP_([[:alnum:]]*).txt','\\1',files)

gl <- grid.layout(nrow=1, ncol=2,widths=1,height=1)
pushViewport(viewport(layout=gl))

vp.TP <- viewport(layout.pos.col=1,layout.pos.row=1)
vp.FP <- viewport(layout.pos.col=2,layout.pos.row=1)
vp <- list("TP"=vp.TP,"FP"=vp.FP)

category <- names(data)
col <- c(brewer.pal(max(3,length(category)),"YlOrRd"))[1:length(category)]
colStr <- paste("fill=c(\"",paste(paste(col),collapse="\",\""),"\")",sep="")

## reading data
for (exp in c("TP","FP"))
{
    pushViewport(vp[[exp]])
    ## get numbers of all intersections
    plotData <- lapply(names(data), function(x) data[[x]][[exp]])
    names(plotData) <- names(data)
    combs <- unlist(lapply(1:length(plotData), function(j) combn(names(plotData), j, simplify = FALSE)),recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    elements <- lapply(combs, function(i) Reduce(intersect,plotData[i]))
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
    grid.text(exp,y=unit(0.05,"npc"),gp=gpar(fontsize=30))
    popViewport()
}
outFile <- file.path(dir,"plots/results/methods_TP_FP.pdf")
dev.print(device=cairo_pdf,outFile,width=40,height=20)