## just load and save data

load("data/mam15.RData")
save(mam15.mt,mam15.ass,mam15.relltest,mam15.aux,
     mam26.mt,mam26.ass,mam26.aux,
     mam105.mt,mam105.ass,mam105.relltest,mam105.aux,file="data/mam15-new.RData")


load("data/lung73.RData")
save(lung73.pvclust,lung73.sb, lung.pvclust,lung.sb,file="data/lung73-new.RData")
