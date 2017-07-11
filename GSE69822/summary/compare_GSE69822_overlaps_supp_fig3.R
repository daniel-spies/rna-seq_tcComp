#!/bin/Rscript
library(VennDiagram)
library(grid)
library(gridBase)
library(RColorBrewer)

path <- "~/Documents/phd/data/simulation/study/results/GSE69822"
comp <-list.dirs(path,recursive=F,full.names=F)
comp <- comp[-grep("shuffle",comp)]
types <- vector("list",length(comp))
names(types) <- comp

## reading in data
for (exp in comp)
{
    files <- grep("edgeR",list.files(file.path(path,exp),recursive=T),invert=T,value=T)
    results <- lapply(files, function(file) read.table(file.path(path,exp,file))[,1])
    names(results) <- sub(".*_([[:alnum:]]*).txt$","\\1",files)
    types[[exp]] <- results
}

# combine pten and pik3ca lists and create new list with wt
combined <- sapply(names(types[[1]]), function(method) unique(c(types[["pten"]][[method]],types[["pik3ca"]][[method]])),simplify=F)
overlaps <- sapply(names(combined), function(meth)list("wt"=types[["wt"]][[meth]],"pten/pik3ca"=combined[[meth]]),simplify=F)

# intersect pertubations
pert_wt <- sapply(names(overlaps),function(x) Reduce(intersect,overlaps[[x]]))

# put all plots in a plotlist and set layout
plotList <- list(overlaps[["DESeq2"]],
                    overlaps[["splineTC"]],
                    overlaps[["maSigPro"]],
                    overlaps[["impulseDE2"]],
                    combined,
                    types[["wt"]],
                    pert_wt)
grid.newpage()
gl <- grid.layout(nrow=6,ncol=12)
pushViewport(viewport(layout=gl))
vp.deseq2       <- viewport(layout.pos.col=1:3,layout.pos.row=1:3)
vp.splineTC     <- viewport(layout.pos.col=4:6,layout.pos.row=1:3)
vp.maSigPro     <- viewport(layout.pos.col=7:9,layout.pos.row=1:3)
vp.impulseDE2   <-viewport(layout.pos.col=10:12,layout.pos.row=1:3)
vp.combined     <- viewport(layout.pos.col=1:4,layout.pos.row=4:6)
vp.wt           <- viewport(layout.pos.col=5:8,layout.pos.row=4:6)
vp.overall      <- viewport(layout.pos.col=9:12,layout.pos.row=4:6)
vp <- list(vp.deseq2,vp.splineTC,vp.maSigPro,vp.impulseDE2,vp.combined,vp.wt,vp.overall)
names(vp) = names(plotList) <- c("DESeq2","splineTC","maSigPro","impulseDE2","pten+pik3ca","wt","pten+pik3ca_vs_wt")

for (plot in names(plotList))
{
    data <- plotList[[plot]]
    data <- data[names(sort(sapply(data,length),decreasing=T))] # reorder data for correct plotting
    category <- names(data)
    names(data) <- 1:length(data)
    pushViewport(vp[[plot]])
    
    # overlapping the pten/pik3ca candidates of each approach
    col <- c(brewer.pal(max(3,length(data)),"YlOrRd"))[1:length(data)]
    colStr <- paste("fill=c(\"",paste(paste(col),collapse="\",\""),"\")",sep="")
    
    ## computing overlaps
    combs <- unlist(lapply(1:length(data), function(j) combn(names(data), j, simplify = FALSE)),recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    elements <- lapply(combs, function(i) Reduce(intersect,data[i]))
    n.elements <- sapply(elements, length)
    
    labels <- paste(category,paste("(n=",n.elements[1:length(data)],")",sep=""),sep="\n")
    catStr <- paste("category=c(\"",paste(labels,collapse="\",\""),"\")",sep="")
    areaStr <- paste(unname(n.elements),collapse=",")
    graphStr <- paste("ind=FALSE,cat.cex=rep(1.5,times=",length(data),"),cat.dist=rep(0.1,times=",length(data),"),cex=rep(1.25,times=",length(elements),"),margin=0.1",sep="")
    param <- paste(areaStr,catStr,colStr,graphStr,sep=",")
    draw <- switch(length(data),"single","pairwise","triple","quad","quintuple")
    grid.draw(eval(parse(text=paste0("draw.",draw,".venn(",param,")"))))
    grid.text(paste(which(names(plotList) == plot),") ",plot," candidates"),y=unit(0.05,"npc"),gp=gpar(fontsize=20))
    popViewport()
    
    # write overlap to list
    overlap <- vector("list",length(data))
    names(overlap) <- length(data):1
    listNames <- category
    symbol <- c()
    
    for (cat in rev(names(elements)))
    {
        cand <- elements[[cat]]
        cand <- cand[!(cand %in% symbol)]
        symbol <- c(symbol,cand)
        naming <- paste(listNames[as.numeric(strsplit(cat,"")[[1]])],collapse="-")

        overlap[[as.character(nchar(cat))]]$symbol <- c(overlap[[as.character(nchar(cat))]]$symbol,cand)
        overlap[[as.character(nchar(cat))]]$count <- c(overlap[[as.character(nchar(cat))]]$count,rep(nchar(cat),length(cand)))
        overlap[[as.character(nchar(cat))]]$method <- c(overlap[[as.character(nchar(cat))]]$method,rep(naming,length(cand)))
    }
    outList <- Reduce(function(u,v) rbind(u,as.data.frame(overlap[[v]])),names(overlap),init=c())
    outFile <- paste0("~/Documents/phd/data/simulation/study/summary/GSE69822/intersection/",plot,"_genes_intersection_list_new.txt")
    write.table(outList,outFile,sep="\t",quote=FALSE,row.names=F,col.names=FALSE)
}
grid.text("Overlap of RNA-seq TC tools on GSE69822 data",gp=gpar(fontsize=40),y=unit(0.98,"npc"))
dev.print(device=cairo_pdf,"~/Documents/phd/data/simulation/study/plots/results/GSE69822/PTEN_PIK3CA_venn_new.pdf",width=60,height=40)