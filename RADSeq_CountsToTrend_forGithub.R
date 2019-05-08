#Webster et al., 2019, RAD-seq starvation resistance manuscript
#This is the script to go from SNP count data to trait data used for GWAS
#The majority of this script was written by Brad Moore

#set working directory so that df.unique.csv and snps.meta.data.csv are in that directory
#read in df.unique.csv and snps.meta.data.csv
#df.unique.csv is the count data, and snps.meta.data.csv describes the replicates
#The ans dataframe in this script will have the separate replicates, and the tmp4 dataframe has the averages

getwd()
setwd("/Users/amywebster/Dropbox (Duke Bio_Ea)/Baugh Lab/baugh_lab_users/Amy/NIL_wormsizer/RADSeq_NIL_paper/RADseq_AmyAnalysis/")

df.unique<-read.csv("df.unique.csv",header = T)
snps.meta.data<-read.csv("snps.meta.data.csv",header = T)

#for the analysis for the RAD-seq paper, we only used indirect replicates 1 and 2
df.unique2<-subset(df.unique,exp.id=="R2D14ind"|exp.id=="R2D1ind"|exp.id=="R2D21ind"|exp.id=="R2D24ind"|exp.id=="R2D7ind"|exp.id=="RID16ind"|exp.id=="RID21ind")
head(df.unique2)

quick.lm <- function(x) {
  ans <- c()
  for (j in unique(x$strain)) {
    for (k in unique(x$exp.id)) {
      tmp <- coef(lm(count ~ total + 0, data=x[x$strain == j & x$exp.id == k,]))
      ans <- rbind(ans, data.frame(slope=tmp["total"], strain=j, exp.id=k))
    }
  }
  ans
}

my.lm <- quick.lm(df.unique2)
my.lm <- merge(my.lm, snps.meta.data, by.x='exp.id', by.y='id')
head(my.lm)

ans <- c()
for (exp in 'indirect') { # look at only indirect results
  tmp <- with(my.lm, my.lm[experiment == exp,])
  for (repl in unique(tmp$replicate.id)) {
    tmp2 <- with(tmp, tmp[tmp$replicate.id == repl,])
    for (s in unique(tmp$strain)) {
      tmp.lm <- lm(slope ~ day, data=with(tmp2, tmp2[strain == s,]))
      ans <- rbind(ans, data.frame(experiment=exp, replicate.id=repl, strain=s, trend=coef(tmp.lm)["day"], r.squared=summary(tmp.lm)$r.squared, rank=0))
    }
    idx <- with(ans, experiment == exp & replicate.id == repl)
    ans$rank[idx] <- (ecdf(ans$trend[idx]))(ans$trend[idx])
  }
}

#Once you have ans, you can use aggregate to get the mean of indirect replicates 1 and 2
trend.lm <- ans
trend.aggr <- aggregate(cbind(trend, r.squared) ~ strain, data=trend.lm, FUN=mean)
trend.lm$strain <- factor(trend.lm$strain, levels=with(trend.aggr, strain[order(trend, decreasing=FALSE)]))
trend.aggr$strain <- factor(trend.aggr$strain, levels=levels(trend.lm$strain))
trend.lm

#df.unique2 <- merge(df.unique, snps.meta.data, by.x='exp.id', by.y='id')

my.r = trend.lm # added 8/29/2018
my.r.aggr <- aggregate(cbind(rank, trend, r.squared) ~ strain, data=my.r, FUN=mean)

#combine trend.r.aggr and my.r.aggr to get the labels in order
tmp3 <- merge(trend.aggr, my.r.aggr, by='strain')
head(tmp3)
with(tmp3, cor(trend.x, trend.y))

names(tmp3) <- c('strain', 'mean.trend.unique', 'mean.squared.unique', 'rank', 'mean.trend.quadOLS', 'mean.rsquared.quadOLS')
#write.csv(tmp3, file='20160408-indirect-results.csv',row.names=FALSE)
tmp4 <- merge(tmp3, with(my.r, my.r[replicate.id == 'R1', c('strain', 'rank')]), by='strain', suffixes=c('.mean', '.R1'))
tmp4 <- merge(tmp4, with(my.r, my.r[replicate.id == 'R2', c('strain', 'rank')]), by='strain', suffixes=c('.mean', '.R2'))
names(tmp4)[8] <- 'rank.R2'
tmp4 <- merge(tmp4, with(trend.lm, trend.lm[replicate.id == 'R1', c('strain', 'trend')], by='strain'))
head(tmp4)

write.csv(tmp4,file='2019-RADSeq-results.csv',row.names = FALSE)

