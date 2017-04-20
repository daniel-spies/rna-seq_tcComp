#!/bin/Rscript
################################### libraries ################################### 
library(VennDiagram)
library(RColorBrewer)
library(grid)
################################### read data ################################### 
path <- "/Users/dspies/Documents/phd/data/simulation/study"
dir <- "results/GSE69822/WT"
wt_tc <- scan("results/GSE69822/WT_TC/results.txt",what="")
inFiles <- list.files(dir,recursive=T)

# replacing wt_tools EBSeq for the paper solution and compare 
files <- list("wt_tools"=inFiles,"wt_tc"=inFiles)
files[["wt_tc"]][1] <- "30M_3rep_6TP_GSE69822_2stepDESeq2.txt"
wt_tools <- lapply(files[["wt_tools"]], function(file) read.table(file.path(dir,file))[,1])
comp <- list("wt_tools"=wt_tools,"wt_tc"=wt_tools)
comp[["wt_tc"]][[1]] <- wt_tc

for (type in names(comp))
{
    dev.new()
    data <- comp[[type]]
    names(data) <- 1:length(data)
    category <- gsub('.*_(.*).txt','\\1',files[[type]])
    col <- c(brewer.pal(max(3,length(category)),"YlOrRd"))[1:length(category)]
    colStr <- paste("fill=c(\"",paste(paste(col),collapse="\",\""),"\")",sep="")

    combs <- unlist(lapply(1:length(data), function(j) combn(names(data), j, simplify = FALSE)),recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    elements <- lapply(combs, function(i) Reduce(intersect,data[i]))
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
    plotDir <- file.path(path,"plots/results/GSE69822")
    system(paste0("mkdir -p ",plotDir))
    dev.print(cairo_pdf,paste0(plotDir,"/",type,"_venn.pdf"),width=20,height=20)
}