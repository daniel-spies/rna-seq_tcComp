#!/bin/Rscript
###################################################### library loading ######################################################
library(biomaRt)
library(FGNet)
################################### functions ################################### 
getHGNC <- function(candidates) {
    symbol <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),filters="ensembl_gene_id",values=candidates,mart=ensembl)
    symbol <- symbol[-which(symbol$hgnc_symbol == "" | duplicated(symbol$ensembl_gene_id)),]
    idx <- candidates %in% symbol$ensembl_gene_id
    candidates[idx] <- symbol$hgnc_symbol

    dIdx <- duplicated(candidates)
    candidates[dIdx] <- renameDuplicates(candidates[dIdx])
    return(candidates)
}
renameDuplicates <- function(nameVector,run=NULL) {
    if(is.null(run)) 
    {
        run <- 2
        nameVector <- paste(nameVector,"_(",run,")",sep="")
    }
    idx <- duplicated(nameVector)
    if(sum(idx) > 0) nameVector[idx] <- renameDuplicates(sub("\\([0-9]*\\)",paste("(",run+1,")",sep=""),nameVector[idx]),run+1)
    else nameVector <- sub("\\([0-9]*\\)",paste("(",run,")",sep=""),nameVector)
    return(nameVector)
}
###################################################### paths ######################################################
dir <- "~/Documents/phd/data/simulation/study/summary/GSE69822"
scriptDir <- "~/Documents/phd/programming/scripts"
outPath <- file.path(dir,"FGNet")
system(paste0("mkdir -p ",outPath))
###################################################### reading data ######################################################
wt <- read.table(file.path(dir,"intersection/wt_genes_intersection_list.txt"))
pten_pik <- read.table(file.path(dir,"intersection/pten+pik3ca_genes_intersection_list.txt"))
wt_t0_pik <- scan("~/Documents/phd/data/simulation/study/results/GSE69822/PIK3A_WT_T0.txt",what="")
wt_t0_pten <- scan("~/Documents/phd/data/simulation/study/results/GSE69822/PTEN_WT_T0.txt",what="")
###################################################### creating lists ######################################################
result_list <- list("DESeq2"=NULL, "splineTC"=NULL, "pairwise"=NULL)
for (type in names(result_list)[-4])
{
    wt_tc <- unique(Reduce(intersect,lapply(unlist(strsplit(type,"-")), function(pat) wt[grep(pat,wt$V3),1])))
    pten_pik_tc <- unique(Reduce(intersect,lapply(unlist(strsplit(type,"-")), function(pat) pten_pik[grep(pat,pten_pik$V3),1])))
    result_list[[type]] <- list(
        "A" = intersect(intersect(wt_tc,pten_pik_tc),intersect(wt_t0_pten,wt_t0_pik)),
        "B" = setdiff(intersect(pten_pik_tc,intersect(wt_t0_pten,wt_t0_pik)),wt_tc),
        "C" = setdiff(intersect(intersect(wt_tc,pten_pik_tc),wt_t0_pik),wt_t0_pten),
        "D" = setdiff(setdiff(intersect(pten_pik_tc,wt_t0_pik),wt_tc),wt_t0_pten),
        "E" = setdiff(intersect(intersect(wt_tc,pten_pik_tc),wt_t0_pten),wt_t0_pik),
        "F" = setdiff(setdiff(intersect(pten_pik_tc,wt_t0_pten),wt_tc),wt_t0_pik),
        "G" = setdiff(setdiff(intersect(wt_tc,pten_pik_tc),wt_t0_pten),wt_t0_pik),
        "H" = setdiff(setdiff(setdiff(pten_pik_tc,wt_tc),wt_t0_pten),wt_t0_pik)
    )
}
combined <- lapply(names(result_list[["DESeq2"]]), function(nList) intersect(result_list[["splineTC"]][[nList]],result_list[["pairwise"]][[nList]]))
names(combined) <- names(result_list[["DESeq2"]])
result_list[["combined_edgeR"]] <- combined
outTable <- rbind(sapply(names(result_list), function(name) sapply(result_list[[name]],length)))

combined <- lapply(names(result_list[["DESeq2"]]), function(nList) intersect(result_list[["splineTC"]][[nList]],result_list[["DESeq2"]][[nList]]))
names(combined) <- names(result_list[["DESeq2"]])
result_list[["combined_DESeq2"]] <- combined
outTable <- rbind(sapply(names(result_list), function(name) sapply(result_list[[name]],length)))

write.table(outTable,file.path(dir,"intersection_table.txt"))
###################################################### FGNet ######################################################
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

for (comp in names(result_list))
{
    path <- file.path(outPath,comp)
    system(paste0("mkdir -p ",path))
    setwd(path)
    
    for (class in names(result_list[[comp]]))
    {
        print(paste0("processing ",comp,":",class))
        data <- result_list[[comp]][[class]]
        label <- getHGNC(data)
        names(label) <- data
            
        tName <- paste("topGO",class,sep="_")
        res_topGO <- fea_topGO(data,geneIdType="ENSEMBL",geneLabels=label,organism="Hs",pValThr = 0.05,jobName=tName)
        outData <- res_topGO$geneTermSets[,c(2,6)]
        outData$Terms <- sub("(GO:[[:digit:]]*):.*","\\1",outData$Terms)
        outFile <- file.path(path,tName,paste(comp,tName,"results.txt",sep="_"))
        write.table(outData,outFile,col.names=F,row.names=F,quote=F,sep="\t")
        system(paste0("python ",scriptDir,"/Revigo.py -i ",outFile," -c 0.50 -o ",tName))
    }
}