#!/bin/Rscript
topN <- 20
path <- "~/Documents/phd/data/simulation/study/summary/GSE69822/FGNet_new"
meth <- list.dirs(path,recursive=F,full.names=F)

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
        topCand$class <- rep(sub("topGO_","",clust),topN)
        # remove enriched categories without sign. genes
        result_clust[[sub("topGO_","",clust)]] <- topCand[complete.cases(topCand),]
    }
    outData <- Reduce(rbind,result_clust)
    outFile <- file.path(path,name,paste0("top",topN,"_GO_summary.csv"))
    write.csv(outData,outFile)
}