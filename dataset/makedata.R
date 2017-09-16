## makedata.R

nb.rell = 100000
nb.pvclust = 10000

date()
library(parallel)
library(pvclust)
library(scaleboot)

length(cl <- makeCluster(detectCores()))

date()

mam15.mt <- read.mt("mam15-files/mam15.mt")
mam15.ass <- read.ass("mam15-files/mam15.ass")
sa <- 9^seq(-1,1,length=13)
mam15.relltest <- relltest(mam15.mt,nb=nb.rell,sa=sa,ass=mam15.ass,cluster=cl)
save(mam15.mt,mam15.ass,mam15.relltest,file="data/mam15.RData")

date()

mam105.mt <- read.mt("mam105-files/mam105.mt")
mam105.ass <- read.ass("mam105-files/mam105.ass")
sa <- 9^seq(-1,1,length=13)
mam105.relltest <- relltest(mam105.mt,nb=nb.rell,sa=sa,ass=mam105.ass,cluster=cl)
save(mam105.mt,mam105.ass,mam105.relltest,file="data/mam105.RData")

date()

data(lung)
sa <- 9^seq(-1,1,length=13) # wider range of scales than pvclust default
lung73.pvclust <- pvclust(lung,r=1/sa,nboot=nb.pvclust,parallel=cl) 
lung73.sb <- sbfit(lung73.pvclust,cluster=cl) # model fitting
save(lung73.pvclust,lung73.sb,file="data/lung73.RData")


date()
quit(save="no")
