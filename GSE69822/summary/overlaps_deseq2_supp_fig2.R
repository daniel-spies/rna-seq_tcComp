#!/bin/Rscript

library(VennDiagram)
library(grid)
library(gridBase)

dev.new()
gl <- grid.layout(nrow=1,ncol=5)
pushViewport(viewport(layout=gl))
vp.deseq2 <- viewport(layout.pos.col=1)
vp.splineTC <- viewport(layout.pos.col=2)
vp.pairwise <- viewport(layout.pos.col=3)
vp.comb_deseq2 <- viewport(layout.pos.col=4)
vp.comb_pairwise <- viewport(layout.pos.col=5)
vp <- list(vp.deseq2,vp.splineTC,vp.pairwise,vp.comb_deseq2,vp.comb_pairwise)
names(vp) <- c("DESeq2","splineTC","pairwise","combined_DESeq2","combined_edgeR")

path <- "~/Documents/phd/data/simulation/study/results/GSE69822"
comp <-list.dirs(path,recursive=F,full.names=F)
types <- vector("list",length(comp))
names(types) <- comp

pten_t0 <- scan(file.path(path,"PTEN_WT_T0.txt"),what="")
pik3ca_t0 <- scan(file.path(path,"PIK3CA_WT_T0.txt"),what="")

## reading in data
for (exp in comp)
{
    files <- list.files(file.path(path,exp),recursive=T)
    results <- lapply(files, function(file) read.table(file.path(path,exp,file))[,1])
    names(results) <- sub(".*_([[:alnum:]]*).txt$","\\1",files)
    types[[exp]] <- results
}

# combine pten and pik3ca lists
combined <- lapply(names(types[[1]]), function(method) unique(c(types[["pten"]][[method]],types[["pik3ca"]][[method]])))
names(combined) <- names(types[[1]])

# combine pten/pik3ca wt overlaps
overlaps <- lapply(names(combined), function(meth){ l <- list("wt"=types[["wt"]][[meth]],"pten/pik3ca"=combined[[meth]]); names(l) <- c("wt",meth); return(l)})
names(overlaps) <- names(combined)

comb_approaches <- lapply(c("DESeq2","pairwise"), function(tool) 
                        lapply(seq_along(overlaps[[tool]]), function(class)
                            intersect(overlaps[[tool]][[class]], overlaps[["splineTC"]][[class]])))

overlaps <- c(overlaps,list("combined_DESeq2"=comb_approaches[[1]]),list("combined_edgeR"=comb_approaches[[2]]))
overlaps <- overlaps[c(1,3,2,4,5)] ## reorder to have the same order as in supp fig. 4

# graphical settings for headings and colours
category <- c("WT EGF-simulated TC","WT vs PTEN KO/PIK3CA","PTEN KO at T=0","PIK3CA at T=0")
areaCol <- c("#B9B9B9","#FFFFFF","#B9B9B9","#B9B9B9")
areaColStr <- paste("fill=c(\"",paste(paste(areaCol),collapse="\",\""),"\")",sep="")
lineCol <- c("#F68B1E","#4573B9","#118340","#010101")
lineColStr <- paste("col=c(\"",paste(paste(lineCol),collapse="\",\""),"\")",sep="")
catCol <- c("#F68B1E","#010101","#4573B9","#118340")
catColStr <- paste("cat.col=c(\"",paste(paste(catCol),collapse="\",\""),"\")",sep="")

for (plot in names(overlaps))
{
    pushViewport(vp[[plot]])
    data <- overlaps[[plot]]
    data <- c(overlaps[[plot]],list(pten_t0),list(pik3ca_t0))
    names(data) <- 1:length(data)
    

    ## computing overlaps
    combs <- unlist(lapply(1:length(data), function(j) combn(names(data), j, simplify = FALSE)),recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    elements <- lapply(combs, function(i) Reduce(intersect,data[i]))
    n.elements <- sapply(elements, length)

    labels <- paste(category,paste("(n=",n.elements[1:length(data)],")",sep=""),sep="\n")
    catStr <- paste("category=c(\"",paste(labels,collapse="\",\""),"\")",sep="")
    areaStr <- paste(unname(n.elements),collapse=",")
    textStr <- paste0("fontfamily=\"Arial\",cat.fontfamily=\"Arial\"")
    graphStr <- paste("ind=FALSE,
                    cat.cex=rep(1.5,times=",length(data),"),
                    cat.dist=rep(0.1,times=",length(data),"),
                    alpha=rep(1,",length(data),"),
                    lwd=rep(5,",length(data),"),
                    cex=rep(1.25,times=",length(elements),"),
                    margin=0.1",sep="")
    param <- paste(areaStr,catStr,areaColStr,lineColStr,catColStr,graphStr,textStr,sep=",")
    draw <- switch(length(data),"single","pairwise","triple","quad","quintuple")
    grid.draw(eval(parse(text=paste0("draw.",draw,".venn(",param,")"))))
    grid.text(plot,y=unit(0.05,"npc"),gp=gpar(fontsize=20))
    popViewport()
}
dev.print(device=cairo_pdf,"~/Documents/phd/data/simulation/study/summary/GSE69822/overlaps_all_venn.pdf",width=100,height=20)