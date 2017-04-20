######################## data input ########################

dir <- "~/Documents/phd/data/simulation/study/countTables"
contr <- read.table(paste(dir,"standard/30M_3rep_4TP_contr.txt",sep="/"))
treat <- read.table(paste(dir,"standard/30M_3rep_4TP_treat.txt",sep="/"))

######################## create samples with lower / higher libSize ########################

libSize <- c (10,20,30,50)
factor_contr <- libSize / round(mean(colSums(contr))/1e6)
factor_treat <- libSize / round(mean(colSums(treat))/1e6)

outData <- list()
for (idx in 1:length(libSize))
{
    repCtrl <- round(contr * factor_contr[idx])
    repTreat <- round(treat * factor_treat[idx])
    outData[[idx]] <- list(contr=repCtrl,treat=repTreat)
}
names(outData) <- paste(libSize,"M_3rep_4TP",sep="")
for (type in names(outData)) for (file in names(outData[[type]])) write.table(outData[[type]][[file]],paste(dir,"/libSize/",type,"_",file,".txt",sep=""))