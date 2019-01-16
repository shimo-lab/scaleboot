## makedata.R

# for making data
nb.rell = 100000
nb.pvclust = 10000
# for testing program
#nb.rell = 1000
#nb.pvclust = 100

date()
library(parallel)
library(pvclust)
library(scaleboot)

length(cl <- makeCluster(detectCores()))

date()

### relltest for mam15
mam15.mt <- read.mt("mam15-files/mam15.mt")
mam15.ass <- read.ass("mam15-files/mam15.ass")
sa <- 9^seq(-1,1,length=13)
mam15.relltest <- relltest(mam15.mt,nb=nb.rell,sa=sa,ass=mam15.ass,cluster=cl)

### aux info for mam15
a <- read.table("mam15-files/mam15.tpl",skip=2)
mam15.tpl <- as.character(a[,1]) # topology
#names(mam15.tpl) <- as.character(a[,2]) # name (actually, the same name as below)
names(mam15.tpl) <- attr(mam15.ass,"trees") # using the names assigned by read.ass
mam15.cld <- as.character(read.table("mam15-files/mam15.cld",skip=4)[,2])
names(mam15.cld) <- attr(mam15.ass,"edges")
mam15.tax <- as.character(read.table("mam15-files/mam15.tax",skip=2)[,2])
mam15.aux <- list(tpl=mam15.tpl, cld=mam15.cld, tax=mam15.tax)

### mam26 ( = mam15 + mam11 ) without relltest
a <- read.table("mam15-files/mam26.tpl",skip=2)
mam26.tpl <- as.character(a[,1]) # topology
names(mam26.tpl) <- as.character(a[,2]) # name
mam26.mt <- read.mt("mam15-files/mam26.mt")
colnames(mam26.mt) <- names(mam26.tpl)
mam26.ass <- read.ass("mam15-files/mam26.ass")
names(mam26.ass)[1:26] <- names(mam26.tpl)
attr(mam26.ass,"trees")[1:26 ] <- names(mam26.tpl)
mam26.aux <- list(tpl=mam26.tpl, cld=mam15.cld, tax=mam15.tax)

date()

### relltest for mam105
mam105.mt <- read.mt("mam15-files/mam105.mt")
mam105.ass <- read.ass("mam15-files/mam105.ass")
sa <- 9^seq(-1,1,length=13)
mam105.relltest <- relltest(mam105.mt,nb=nb.rell,sa=sa,ass=mam105.ass,cluster=cl)


### aux info for mam105
a <- read.table("mam15-files/mam105.tpl",skip=2)
mam105.tpl <- as.character(a[,1]) # topology
#names(mam15.tpl) <- as.character(a[,2]) # name (actually, the same name as below)
names(mam105.tpl) <- attr(mam105.ass,"trees") # using the names assigned by read.ass
mam105.cld <- as.character(read.table("mam15-files/mam105.cld",skip=4)[,2])
names(mam105.cld) <- attr(mam105.ass,"edges")
mam105.tax <- as.character(read.table("mam15-files/mam105.tax",skip=2)[,2])
mam105.aux <- list(tpl=mam105.tpl, cld=mam105.cld, tax=mam105.tax)

### save the result for mam15 and mam11
save(mam15.mt,mam15.ass,mam15.relltest,mam15.aux,
     mam26.mt,mam26.ass,mam26.aux,
     mam105.mt,mam105.ass,mam105.relltest,mam105.aux,file="data/mam15.RData")

date()

data(lung)
sa <- 9^seq(-1,1,length=13) # wider range of scales than pvclust default
lung73.pvclust <- pvclust(lung,r=1/sa,nboot=nb.pvclust,parallel=cl) 
lung73.sb <- sbfit(lung73.pvclust,cluster=cl) # model fitting
save(lung73.pvclust,lung73.sb,file="data/lung73.RData")


date()
quit(save="no")
