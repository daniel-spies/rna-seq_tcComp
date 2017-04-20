#!/bin/Rscript

library(VennDiagram)
library(grid)
library(gridBase)

path <- "~/Documents/phd/data/simulation/study/results/GSE69822"
comp <- grep("shuffle",list.dirs(path,recursive=F,full.names=F),invert=T,value=T)
types <- vector("list",length(comp))
names(types) <- comp

pten_t0 <- scan(file.path(path,"PTEN_WT_T0.txt"),what="")
pik3ca_t0 <- scan(file.path(path,"PIK3CA_WT_T0.txt"),what="")

## reading in data
for (exp in comp)
{
    files <- grep("edgeR",list.files(file.path(path,exp),recursive=T),value=T,invert=T)
    results <- lapply(files, function(file) read.table(file.path(path,exp,file))[,1])
    names(results) <- sub(".*_([[:alnum:]]*).txt$","\\1",files)
    types[[exp]] <- results
}

# combine pten and pik3ca lists and create new list with wt
combined <- sapply(names(types[[1]]), function(method) unique(c(types[["pten"]][[method]],types[["pik3ca"]][[method]])),simplify=F)
overlaps <- sapply(names(combined), function(meth)list("wt"=types[["wt"]][[meth]],"pten/pik3ca"=combined[[meth]]),simplify=F)

# intersect candidates of each type and add to plotList
comb_approaches <- sapply(names(types), function(exp) Reduce(intersect,types[[exp]]))
comb_list <- list("wt" = comb_approaches[["wt"]], "pten/pik3ca" = unique(c(comb_approaches[["pten"]],comb_approaches[["pik3ca"]])))
overlaps <- c(overlaps,"combined"=list(comb_list))

# graphical settings for headings and colours
category <- c("WT EGF-simulated TC","WT vs PTEN KO/PIK3CA","PTEN KO at T=0","PIK3CA at T=0")
areaCol <- c("#B9B9B9","#FFFFFF","#B9B9B9","#B9B9B9")
areaColStr <- paste("fill=c(\"",paste(paste(areaCol),collapse="\",\""),"\")",sep="")
lineCol <- c("#F68B1E","#4573B9","#118340","#010101")
lineColStr <- paste("col=c(\"",paste(paste(lineCol),collapse="\",\""),"\")",sep="")
catCol <- c("#F68B1E","#010101","#4573B9","#118340")
catColStr <- paste("cat.col=c(\"",paste(paste(catCol),collapse="\",\""),"\")",sep="")

# graphical settings for layout
dev.new()
gl <- grid.layout(nrow=2,ncol=2)
pushViewport(viewport(layout=gl))
vp.deseq2 <- viewport(layout.pos.col=1,layout.pos.row=1)
vp.splineTC <- viewport(layout.pos.col=2,layout.pos.row=1)
vp.maSigPro<- viewport(layout.pos.col=1,layout.pos.row=2)
vp.combined <- viewport(layout.pos.col=2,layout.pos.row=2)
vp <- list(vp.deseq2,vp.splineTC,vp.maSigPro,vp.combined)
names(vp) <- c("DESeq2","splineTC","maSigPro","combined")

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
    textStr <- paste0("fontfamily=\"Arial\",cat.fontfamily=\"Helvetica\"")
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
dev.print(device=cairo_pdf,"~/Documents/phd/data/simulation/study/plots/results/GSE69822/overlaps_all_venn_new.pdf",width=20,height=20)