library(ImpulseDE2)

###################################################### data input ######################################################
args <- commandArgs(trailingOnly = TRUE)if(length(args) < 2) args <- c("--help")
if("--help" %in% args) 
{
    cat("Arguments:\n--args inFile1 inFile2 inFileN outDir\n\nFilenaming contains parameters: e.g. xM_yrep_zTP_[contr/treat].txt")
    q(save="no")
}

outPath <- args[length(args)]
system(paste0("mkdir -p ",outPath))
for (idx in 1:(length(args)-1))
{
    # parse parameters
    print(paste("analyzing:",basename(args[idx])))
    par <- strsplit(basename(args[idx]),"_")[[1]]
    repl <- as.numeric(strsplit(par[2],"rep")[[1]])
    TP <- as.numeric(strsplit(par[3],"TP")[[1]])
    time <- ((1:TP)-1)*3
    pVal <- 0.01
    header <- paste(rep(((1:TP)-1)*3,each=repl),"h-",rep(1:repl),sep="")

    # read data
    contr <- read.table(args[idx],header=T,stringsAsFactors=F)
    treat <- read.table(sub("contr","treat",args[idx]),header=T,stringsAsFactors=F)
    rownames(treat) <- rownames(contr)

    # merge data into 1 data.framea
    inData <- cbind(contr,treat)
    type <- c("control","case")
    colnames(inData) <- paste(rep(type,each=length(header)),rep(header,2),sep="_")
    inData <- as.matrix(inData[rowSums(inData)>0,]) # otherwise crash
    
    # specify experimental design
    design <- data.frame("Sample"=colnames(inData),
                         "Condition"=rep(type,each=TP*repl),
                         "Time"=rep(rep(time,each=repl),2),
                         "Batch"=rep("B_NULL",ncol(inData)), 
                         row.names=colnames(inData))
    # DEG analysis
    impulse_results <- runImpulseDE2(matCountData = inData,
                                    dfAnnotation =design,
                                    boolCaseCtrl = TRUE,
                                    scaNProc = 20,
                                    scaQThres = pVal,
                                    boolIdentifyTransients = TRUE)

    outFile <- sub("(contr|treat).txt","ImpulseDE2.txt",basename(args[idx]))
    write.table(impulse_results$dfImpulseDE2Results[,c(1,3)],file.path(outPath,outFile),row.names=F,col.names=F,quote=F,sep="\t")
}