#!/bin/Rscript
# comparing top 20 GO terms of each method and intersection
library(gridExtra)
library(ggplot2)

topN <- 10
path <- "~/Documents/phd/data/simulation/study/summary/GSE69822/FGNet_new2"
meth <- list.dirs(path,recursive=F,full.names=F)

result_list <- vector("list",length(meth))
names(result_list) <- meth

# extract information
for (name in meth)
{
    clusters <- list.dirs(file.path(path,name),recursive=F,full.names=F)
    result_clust <- vector("list",length(clusters))
    names(result_clust) <- sub("topGO_","",clusters)
    for (clust in clusters)
    {
        #get number of significant genes in each cluster
        nSig <- read.table(file.path(path,name,clust,paste0(clust,".txt")),sep="\t",header=T)[,c(2,4)]
        rownames(nSig) <- sub("(GO:[[:digit:]]*):.*","\\1",nSig$Terms)
        
        # sort by p-value and only use collapsed groups
        enrichment <- read.csv(file.path(path,name,clust,"REVIGO.csv"))
        enrichment <- enrichment[order(enrichment$log10.p.value),]
        enrichment <- enrichment[enrichment$eliminated == 0,c("term_ID","description","log10.p.value")]
        topCand <- enrichment[1:topN,]
        topCand$sig <- nSig[topCand$term_ID,"Significant"]
        # remove enriched categories without sign. genes
        result_clust[[sub("topGO_","",clust)]] <- topCand[complete.cases(topCand),]
    }
    result_list[[name]] <- result_clust
}

# plotting of tables
plotList <- list("space"=ggplot() + theme_bw() + theme(plot.margin=unit(c(rep(1,4)),"npc"))) ## start with empty space
plotNames <- c("DESeq2","impulseDE2","splineTC","maSigPro",
               "DESeq2_impulseDE2","DESeq2_splineTC","DESeq2_maSigPro","impulseDE2_splineTC","impulseDE2_maSigPro",
               "DESeq2_splineTC_impulseDE2","maSigPro_splineTC_impulseDE2",
               "all")
colorList <- list("#73BF44","#8D6914","#D66EAA","#EAE60F",
                  "#818F29","#AA927D","#C2D921","#BE6C78","#CBBC11",
                  "#C7781D","#FF3333",
                  "#C21807")
names(colorList) <- plotNames

# create heading
text.size=24
plotList <- c(plotList,sapply(names(result_list[[1]]), function(label) 
    ggplot() + annotate("text", x = 4, y = 25, size=text.size, label = label) + theme_bw() + theme_void(),
    simplify=F))

print(length(plotList))
for (label in plotNames)
{
    plotList[[length(plotList)+1]] <- ggplot() + annotate("text", x = 4, y = 25, size=text.size, label = label) + theme_bw() + theme_void()
    for (clust in names(result_list[[label]]))
    {
        d <- result_list[[label]][[clust]]
        d$description <- factor(d$description,levels=rev(d$description))
        # creating labels for first and last element
        labels <- d
        labels$log10.p.value <- round(d$log10.p.value,2)
        labels$log10.p.value[2:(nrow(labels)-1)] <- NA
        # creating plots
        plotList[[length(plotList)+1]] <- 
            ggplot(d, aes(description, sig, fill = log10.p.value)) +
                geom_bar(stat = "identity") +
                coord_flip() +
                geom_text(data=labels, aes(label=log10.p.value), vjust=0.5, hjust=0) +
                scale_fill_gradient(high = "grey40", low = colorList[[label]]) +
                theme_bw() +
                theme(axis.line = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                     panel.background = element_blank(),
                     legend.position = "none",
                     axis.title.x = element_blank(),axis.title.y = element_blank(),
                     axis.text.y = element_text(size=text.size))
    }
}

# transform ggplots into grobPlots
grobPlots <- lapply(1:length(plotList), function(i) ggplotGrob(plotList[[i]]))
ncol=length(result_list[[1]]) + 1
nrow=c(round(length(grobPlots)/ncol/2), (length(plotNames) + 1) - round(length(grobPlots)/ncol/2))

# split onto two pdf pages
pdf("~/Documents/phd/data/simulation/study/plots/results/GSE69822/GSE69822_GO_results_new2_top10.pdf",width=200,height=200)
grid.arrange(grobs=grobPlots[1:(ncol*nrow[1])],ncol=ncol,nrow=nrow[1],as.table=T)
grid.arrange(grobs=grobPlots[(ncol*nrow[1]+1):length(plotList)],ncol=ncol,nrow=nrow[2],as.table=T)
dev.off()