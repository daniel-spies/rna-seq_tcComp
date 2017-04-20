#!/bin/Rscript

library(VennDiagram)
library(grid)
library(RColorBrewer)

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

# plotting Venn diagram for DESe2 used for figure 3
plot<- "DESeq2"
data <- c(overlaps[[plot]],list("PTEN_T0"=pten_t0),list("PIK3CA_T0"=pik3ca_t0))

category <- names(data)
names(data) <- 1:length(data)

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
dev.print(device=cairo_pdf,"~/Documents/phd/data/simulation/study/summary/GSE69822/overlaps_DESeq2_venn.pdf",width=20,height=20)