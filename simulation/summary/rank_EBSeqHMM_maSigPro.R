path <- "/Users/dspies/Documents/phd/data/simulation/study"
###########################################################################
# RENAME "all" files first if you want to keep the original results !!!!! #
######################## maSigPro #########################################
deg <- read.table(file.path(path,"results/maSigPro/30M_3rep_4TP_maSigPro.txt"))
all <- read.table(file.path(path,"results_noCutoff/maSigPro/30M_3rep_4TP_maSigPro_all.txt"))
truth <- read.table(file.path(path,"stats/SIM_truth.txt"),header=T)

# add missing genes to maSigPro all
all <- rbind(all,data.frame(V1=setdiff(truth$feature,all$V1),V2=1,V3=1))

# split groups wether they are DEG or not
notDEG <- all[which(!(all$V1 %in% unlist(deg))),]

# group genes by their PP and sort decreasing
#rankDEG <- aggregate( V1 ~ V2,deg, c)
#rankDEG <- rankDEG[order(rankDEG$V2),]
#rank <- nrow(rankDEG):1
#outData <- Reduce(rbind,sapply(1:length(rank), function(idx)cbind(rankDEG$V1[[idx]],rank[[idx]])))
#outData <- rbind(outData, cbind(notDEG$V1,1))
outData <- rbind(deg[,1:2],cbind(notDEG[,1],0.999))
write.table(outData,file.path(path,"results_noCutoff/maSigPro/30M_3rep_4TP_maSigPro.txt"),row.names=F,col.names=F,quote=F)

######################## EBSeqHMM ########################
deg <- read.table(file.path(path,"results/EBSeqHMM/30M_3rep_4TP_EBSeqHMM.txt"))
all <- read.table(file.path(path,"results_noCutoff/EBSeqHMM/30M_3rep_4TP_EBSeqHMM_all.txt"))

# add missing genes to maSigPro all
all <- rbind(all,data.frame(V1=setdiff(truth$feature,all$V1),V2=0))

# split groups wether they are DEG or not
notDEG <- all[which(!(all$V1 %in% deg$V1)),]

# group genes by their PP and sort decreasing
rankDEG <- aggregate( V1 ~ V3,deg, c)
rankDEG <- rankDEG[order(rankDEG$V3,decreasing=T),]

rankNotDEG <- aggregate( V1 ~ V2,notDEG, c)
rankNotDEG <- rankNotDEG[order(rankNotDEG$V2,decreasing=T),]

#combine lists and cast back into long format
colnames(rankNotDEG) = colnames(rankDEG) <- c("score","gene")
newOrder <- rbind(rankDEG,rankNotDEG)

rank <- nrow(newOrder):1
outData <- Reduce(rbind,sapply(1:length(rank), function(idx)cbind(newOrder$gene[[idx]],rank[[idx]])))
write.table(outData,file.path(path,"results_noCutoff/30M_3rep_4TP_EBSeqHMM.txt"),row.names=F,col.names=F,quote=F)

######################## FunPat ########################
all <- read.table(file.path(path,"results_noCutoff/FunPat/30M_3rep_4TP_FunPat_all.txt"))
all$V1 <- sub("gene_","",all$V1)

# add missing genes to FunPat all
all <- rbind(data.frame(V1=setdiff(truth$feature,all$V1),V2=1),all)
write.table(all,file.path(path,"results_noCutoff/FunPat/30M_3rep_4TP_FunPat.txt"),row.names=F,col.names=F,quote=F)