#!/bin/Rscript
## set paths and load identifier files
dir <- "~/Documents/phd/data/simulation/study"
files <- list.files(file.path(dir,"results"),pattern="30M_3rep_4TP.*txt",recursive=T)
files <- files[-grep("(noise|GSE|timeSeq|nsgp|EBSeqHMM)",files)]
DEG <- scan(file.path(dir,"stats/DEG_IDs.txt"),what="",sep="\n")
NEG <- scan(file.path(dir,"stats/NEG_IDs.txt"),what="",sep="\n")
SIM <- read.table(file.path(dir,"stats/SIM_IDs.txt"),header=F,stringsAsFactors=F)

# creating list of TP / FP
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

## reading data
for (exp in c("TP","FP"))
{
    ## get numbers of all intersections
    overlapData <- lapply(names(data), function(x) data[[x]][[exp]])
    names(overlapData) <- 1:length(data)
    combs <- unlist(lapply(1:length(overlapData), function(j) combn(names(overlapData), j, simplify = FALSE)),recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    elements <- lapply(combs, function(i) Reduce(intersect,overlapData[i]))
    n.elements <- sapply(elements, length)
    
    # get number of genes for each overlap 1:length(data)
    outData <- t(sapply(1:length(data), function(tool) {
                    idx <- grep(as.character(tool),names(elements),value=T)
                    used <- c()
                    outVector <- c()
                    for (i in max(nchar(idx)):1)
                    {
                        overlap_genes <- unique(unlist(elements[idx[which(nchar(idx) == i)]]))
                        overlap_genes <- overlap_genes[!(overlap_genes %in% used)]
                        used <- c(used,overlap_genes)
                        outVector <- c(outVector,length(overlap_genes))
                    }
                    return(outVector)
                }))
    colnames(outData) <- length(data):1
    rownames(outData) <- names(data)
    if (exp == "TP") outData <- outData[do.call(order,args=c(as.data.frame(outData),decreasing=T)),]
        else outData <- outData[do.call(order,args=c(rev(as.data.frame(outData)),decreasing=F)),]
    outFile <- file.path(dir,paste0("stats//methods_overlap_",exp,"_top5.txt"))
    write.table(outData,outFile,row.names=T,col.names=T,quote=F)
}