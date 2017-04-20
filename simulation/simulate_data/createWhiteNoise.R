#!/bin/Rscript
# create noisy data

path <- "~/Documents/phd/data/simulation/study"
contr <- read.table(file.path(path,"run_data/30M_3rep_4TP_contr.txt"),row.names=1)
treat <- read.table(file.path(path,"run_data/30M_3rep_4TP_treat.txt"),row.names=1)
rownames(treat) <- rownames(contr)

noise <- c(0,0.05,0.1,0.15,0.2)
for (w.noise in noise)
{
    n <- length(unlist(contr))
    contr.wn <- contr * rnorm(n,1,w.noise)
    write.table(contr.wn,paste0(path,"/countTables/noise/30M_3rep_4TP_",w.noise,"wn_contr.txt"))
    treat.wn <- treat * rnorm(n,1,w.noise)
    write.table(treat.wn,paste0(path,"/countTables/noise/30M_3rep_4TP_",w.noise,"wn_treat.txt"))
}