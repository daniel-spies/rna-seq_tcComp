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
outPath <- file.path(dir,"FGNet_new2")
system(paste0("mkdir -p ",outPath))
###################################################### reading data ######################################################
wt <- read.table(file.path(dir,"intersection/wt_genes_intersection_list_new.txt"))
pten_pik <- read.table(file.path(dir,"intersection/pten+pik3ca_genes_intersection_list_new.txt"))
wt_t0_pik <- scan("~/Documents/phd/data/simulation/study/results/GSE69822/PIK3CA_WT_T0.txt",what="")
wt_t0_pten <- scan("~/Documents/phd/data/simulation/study/results/GSE69822/PTEN_WT_T0.txt",what="")
###################################################### creating lists ######################################################
result_list <- list("DESeq2"=NULL, "splineTC"=NULL,"maSigPro"=NULL,"impulseDE2"=NULL)
for (type in names(result_list))
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
###### add combined list to results
## DESeq2 + [impulseDE2 / splineTC / maSigPro] methods (pairwise combination)
combined_spline <- sapply(c("impulseDE2","splineTC","maSigPro"), function(pair)
    sapply(LETTERS[1:8], function(nList) 
        intersect(result_list[["DESeq2"]][[nList]],result_list[[pair]][[nList]]),
        simplify=F),
simplify=F)
names(combined_spline) <- paste0("DESeq2_",names(combined_spline))

## impulseDE2 + [ splineTC / maSigPro] methods (TC combination)
combined_spline2 <- sapply(c("splineTC","maSigPro"), function(pair)
    sapply(LETTERS[1:8], function(nList) 
        intersect(result_list[["impulseDE2"]][[nList]],result_list[[pair]][[nList]]),
        simplify=F),
simplify=F)
names(combined_spline2) <- paste0("impulseDE2_",names(combined_spline2))

## [DESeq + 2 best TC tools / all 3 best TC tools] combined
combined_spline_maSigPro <- sapply( c("maSigPro","DESeq2"), function(pair)
    sapply(LETTERS[1:8], function(nList){
        intersectList <- list(result_list[["splineTC"]][[nList]],
                              result_list[["impulseDE2"]][[nList]],
                              result_list[[pair]][[nList]])
        Reduce(intersect,intersectList)
    },simplify=F),
simplify=F)
names(combined_spline_maSigPro) <- paste0(names(combined_spline_maSigPro),"_splineTC_impulseDE2")

#all together
combined_all <- sapply(LETTERS[1:8], function(nList) Reduce(intersect,sapply(result_list,"[[",nList)))

result_list <- c(result_list,combined_spline,combined_spline2,combined_spline_maSigPro,"all"=list(combined_all))
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