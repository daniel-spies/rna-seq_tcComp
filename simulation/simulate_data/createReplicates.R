######################## data input ########################

dir <- "~/Documents/phd/data/simulation/study/countTables"
contr <- read.table(paste(dir,"raw/contr_1.txt",sep="/"))
treat <- read.table(paste(dir,"raw/treat_1.txt",sep="/"))

######################## create samples with less replicates ########################

nRep <- c(2,3,5)
outData <- list()
cols <- c(1,11,21,31)

for (idx in 1:length(nRep))
{
    repCtrl <- contr[,unlist(lapply(cols,function(x) seq(x,x+nRep[idx]-1)))]
    repTreat <- treat[,unlist(lapply(cols,function(x) seq(x,x+nRep[idx]-1)))]
    outData[[idx]] <- list(contr=repCtrl,treat=repTreat)
}
names(outData) <- paste("30M_",nRep,"rep_4TP",sep="")
for (type in names(outData)) for (file in names(outData[[type]])) write.table(outData[[type]][[file]],paste(dir,"/replicates/",type,"_",file,".txt",sep=""))