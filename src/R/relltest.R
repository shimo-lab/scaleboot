##
##  scaleboot: R package for multiscale bootstrap
##  Copyright (C) 2006 Hidetoshi Shimodaira
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
##
######################################################################
### MAIN: RESAMPLING


### scaleboot
## dat: data matrix (array of row vectors)
## nb: vector of numbers of replicates
## sa: vector of sigma^2
## fun(x,i,parm): calculate a replicate from indices
## fun(x,w,parm): calculate a replicate from weights
## parm: parameter to be passed to fun
## count: whether returns only counts or statistics
## weight: whether resample via weights or indecies
## cluster: for parallel computing (snow package)
##
## output: list(nb,sa,nsize,nsizes,cnt,bps,stat)
## sa: recomputed scales
## nb: as used
## bps: matrix of bootstrap probabilities (only for count=T)
## stat: list of replicates, i.e., raw output of fun (only for count=F)

scaleboot <- function(dat,nb,sa,fun,parm=NULL,count=TRUE,weight=TRUE,cluster=NULL,onlyboot=FALSE,seed=NULL,...) {
  ## init
  nsize <- nrow(dat)  # (n) sample size of dataset
  nscale <- length(sa) # number of scales
  nb <- rep(nb,length=nscale)
  if(!is.null(seed)) {
    set.seed(seed)
    if(!is.null(cluster)) {
      seeds <- seed+seq(along=cluster)
      parLapply(cluster, as.list(seeds), set.seed)
    }
  }

  ## recompute sa
  nsizes <- round(nsize/sa)  # (n') sample sizes of replicates
  sa <- nsize/nsizes

  ## resampling
  if(is.null(cluster)) {
    ## for single cpu
    z <- scaleboot.node(dat,nb,nsizes,fun,parm,count,weight)
  } else {
    ## for parallel cpu's ##
    fun2 <- fun;  dat2 <- dat; parm2 <- parm # for exporting
    ncl <- length(cluster)
    nbs <- splitnum(nb,ncl)
    arg <- vector("list",ncl)
    for(i in 1:ncl) arg[[i]] <- nbs[i,]
    zz <- parLapply(cluster,arg,
              function(nb1) scaleboot.node(dat2,nb1,nsizes,fun2,parm2,count,weight))
    ## collect results from cpu's
    z <- vector("list",nscale)
    if(count) { ## for counts
      for(i in 1:nscale) {
        ## summing up the count vectors
        z[[i]] <- 0
        for(j in seq(1,length=ncl)) {
          if(!is.numeric(y <- zz[[j]][[i]])) cat(y)
          else z[[i]] = z[[i]] + y
        }
      }
    } else {  ## for stats
      for(i in 1:nscale) {
        ## results are joined together
        k <- 1
        for(j in 1:ncl) {
          len <- length(zz[[j]][[1]])
          z[[i]][seq(k,length=len)] <- zz[[j]][[1]]
          zz[[j]][[1]] <- NULL  ## to shift left the elemnts of zz[[j]]
          k <- k + len
        }
      }
    }
    ## end of parallel cpu's ##
  }
  ## return value
  if(count) {
    cnt <- asmat(simplist(z))
    bps <- sweep(cnt,2,nb,"/")
    if(nrow(bps)==1) bps <- bps[1,]  # as.vector
    z=NULL
  } else bps=NULL
  if(onlyboot) {
    return(list(stat=z,bps=bps,nb=nb,sa=sa))
  } else if(count) {
    sbfit(bps,nb,sa,cluster=cluster,...)
  } else {
    sbconf(z,sa,cluster=cluster,...)
  }
}

### internal resampling function
## dat: data matrix (array of row vectors)
## nb: vector of numbers of replicates
## nsizes: vector of sample sizes
## fun(x,i,parm): calculate a replicate from indices
## fun(x,w,parm): calculate a replicate from weights
## parm: parameter to be passed to fun
## count: whether returns only counts or statistics
## weight: whether resample via weights or indecies
scaleboot.node <- function(dat,nb,nsizes,fun,parm,count,weight) {
  nsize <- nrow(dat)  # n
  w0 <- rep(1,nsize)
  nscale <- length(nsizes) # number of scales
  nb <- rep(nb,length=nscale)
  z <- vector("list",nscale)
  for(i in seq(length=nscale)) {
    nsize2 <- nsizes[i]  # n'
    y <- vector("list",nb[i])
    if(weight) { ## multinomial distribution
      for(j in seq(length=nb[i])) y[[j]] <- fun(dat,
                     as.vector(rmultinom(1,nsize2,w0)),parm)
    } else { ## resampling indices
      for(j in seq(length=nb[i])) y[[j]] <- fun(dat,
                     sample(nsize,nsize2,replace=T),parm)
    }
    if(count) {
      y <- apply(trmat(simplist(y)),2,sum)
    }
    z[[i]] <- y
  }
  z
}

### examples for "fun" in scaleboot.node
## calculate test statistics for selection problem
## via the RELL resamling method
##
## fun(x,i,parm) or fun(x,w,parm)
## x : data matrix (each row is an item)
## i : indecies to be sampled
## w : weights for items;
##   w=(1,...,1) corresponds to the original data
## parm : parameter to be passed to fun
##        (this formal parameter is needed even if unused)
##
## There are four types of functions determied by
##
## * outputs
##     statistics : count=F
##     counts : count=T
## * selected by
##     indices : weight=F
##     weights : weight=T
##
## The default is (count=T, weight=T), i.e., type4 below.
##
## type1: (count=F, weight=F)
#stati.max <- function(x,i,a) maxdif(sumrow(x,i))
## type2: (count=F, weight=T)
#statw.max <- function(x,w,a) maxdif(wsumrow(x,w))
## type3: (count=T, weight=F)
#counti.max <- function(x,i,a) maxdif(sumrow(x,i)) <= 0
## type4: (count=T, weight=T) this gives the best performance.
#countw.max <- function(x,w,a) maxdif(wsumrow(x,w)) <= 0
## a generalization of "type4" above
## association max
## this is equivalent to countw.max if a=NULL
countw.assmax <- function(x,w,ass) {
  y <- maxdif(wsumrow(x,w)) <= 0 # countw.max
  if(is.null(ass)) y
  else {
    z <- vector("logical",length(ass))
    for(i in seq(along=ass)) z[i] <- any(y[ass[[i]]])
    z
  }
}
## yet other examples for Shimodaira-Hasegawa test
##  obs=observed statistics
countw.shtest <- function(x,w,obs)  maxdif(wsumrow(x,w)) >= obs
## parm=list(ass,obs)
countw.shtestass <- function(x,w,assobs)
  unlist(assmaxdif(wsumrow(x,w),assobs$ass)) >= assobs$obs


######################################################################
### MAIN: TESTING FOR SELECTION PROBLEM via RELL METHOD


##
## extracting elements
##

sbextpv <- function(x,i) {
  list(pv=x$pv[i],pe=x$pe[i],nb=x$nb)
}

"[.relltest" <- function(x, i, ...)
  structure(NextMethod("["),
            stat = attr(x,"stat")[i],
            shtest = sbextpv(attr(x,"shtest"),i)
            )


### FRONT-END
##
## dat: data matrix (prodeuced by read.mt)
## nb: number of replicates
## ass: association list (produced by read.ass)
## cluster: for parallel computing (snow package)

relltest <- function(dat,nb=10000,sa=9^seq(-1,1,length=13),
                     ass=NULL,cluster=NULL,nofit=FALSE,models=NULL,
                     seed=100){

  if(ncol(dat)<2) stop("should be more than two trees")
  if(is.null(ass)) na <- colnames(dat) else na <- names(ass)
  ans <- scaleboot(dat,nb,sa,countw.assmax,ass,cluster=cluster,
                   names.hp=na,nofit=nofit,models=models,seed=seed)
  ans2 <- rell.shtest(dat,nb,ass,cluster,seed=NULL)
  attr(ans,"stat") <- ans2$stat
  attr(ans,"shtest") <- ans2$shtest
  class(ans) <- c("relltest",class(ans))
  ans    
}

table.relltest <- function(x,...) {
  mycatpval <- function(x) catpval(x$pv,x$pe)$value
  y <- cbind(myformat(c(pi,attr(x,"stat")),digits=2)[-1],
             mycatpval(attr(x,"shtest")))
  colnames(y) <- c("stat","shtest")
  y
}

print.relltest <- function(x,...) {
  cat("\nTest Statistic, and Shimodaira-Hasegawa test:\n")
  catmat(table.relltest(x,...))
  NextMethod("print")
}

### Shimodaira-Hasegawa test
##
## dat: data matrix (prodeuced by read.mt)
## nb: number of replicates (scalar)
## ass: association list (produced by read.ass)
## cluster: for parallel computing (snow package)

rell.shtest <- function(dat,nb,ass=NULL,cluster=NULL,seed=NULL) {
  nb <- round(mean(nb))
  lik <- apply(dat,2,sum)
  n <- nrow(dat)
  z <- sweep(dat,2,lik/n) # centering
  if(is.null(ass)) { # for comparing trees (SHi)
    tobs <- maxdif(lik) # observed statistics
    sim <-scaleboot(z,nb,1,countw.shtest,tobs,cluster=cluster,
                    onlyboot=TRUE,seed=seed)
    pv <- as.vector(sim$bps) # p-value
    names(pv) <- names(lik)
  } else { # for comparing trees and edges (SHi and SHe)
    aobs <- assmaxdif(lik,ass)
    tobs <- sapply(aobs,min)
    pa <- list(ass=ass,obs=unlist(aobs))
    sim <- scaleboot(z,nb,1,countw.shtestass,pa,cluster=cluster,
                     onlyboot=TRUE,seed=seed)
    pv0 <- as.vector(sim$bps)
    j2 <- cumsum(sapply(ass,length)); j1 <- 1+c(0,j2[-length(j2)])
    pv1 <- aobs
    for(i in seq(along=pv1)) pv1[[i]][] <- pv0[j1[i]:j2[i]]
    pv <- sapply(pv1,max)
  }
  
  list(stat=tobs,shtest=list(pv=pv,pe=sebp(pv,nb),nb=nb))
}


## x: edge2tree
## output: tree2edge
revmap <- function(x,lab="t") {
  n <- max(unlist(x))
  y <- vector("list",n)
  if(!is.null(lab)) names(y) <- paste(lab,seq(along=y),sep="")
  for(i in seq(along=x)) {
    for(j in x[[i]]) y[[j]] <- c(y[[j]],i)
  }
  y
}

##########################################################################
### PHYLOGENETIC ANALYSIS WITH CONSEL
### 
### prepare table for phylogenetic inference (sort by stat)
## trees, edges: outout from relltest
## edge2tree: ass information from treeass (consel)
## treename: vector of tree names
## edgename: vector of edge names
prepare.table.tree <- function(trees,edges,edge2tree,treename=NULL,edgename=NULL,taxaname=NULL,mt=NULL) {
  tree2edge <- revmap(edge2tree)
  
  stat.tree <- attr(trees,"stat") # likelihood
  stat.edge <- attr(edges,"stat") # likelihood
  
  order.tree <- order(stat.tree)
  names(order.tree) <- paste("T",seq(along=order.tree),sep="")
  order.edge <- order(stat.edge)
  names(order.edge) <- paste("E",seq(along=order.edge),sep="")
  invorder.tree <- invPerm(order.tree)
  names(invorder.tree) <- names(trees)
  invorder.edge <- invPerm(order.edge)
  names(invorder.edge) <- names(edges)
  
  trees <- trees[order.tree]; names(trees) <- names(order.tree)
  edges <- edges[order.edge]; names(edges) <- names(order.edge)
  if(!is.null(treename)) {
    treename <- treename[order.tree]; names(treename) <- names(order.tree)
  }
  if(!is.null(edgename)) {
    edgename <- edgename[order.edge]; names(edgename) <- names(order.edge)
  }
  if(!is.null(mt)) {
    mt <- mt[,order.tree]; colnames(mt) <- names(order.tree)
  }
  
  new.edge2tree <- lapply(edge2tree, function(a) unname(invorder.tree[a]))[order.edge]
  names(new.edge2tree) <- names(order.edge)
  new.tree2edge <- lapply(tree2edge, function(a) unname(invorder.edge[a]))[order.tree]
  names(new.edge2tree) <- names(order.edge)
  
  list(trees=trees,edges=edges,edge2tree=new.edge2tree, tree2edge=new.tree2edge, 
       order.tree=order.tree,order.edge=order.edge,invorder.tree=invorder.tree,invorder.edge=invorder.edge,
       treename=treename, edgename=edgename, taxaname=taxaname,mt=mt)
}

### prepare tables for phylogenetic analysis
### 
## x: output from reorder.tree
table.tree <- function(x,k=2) {
  opt.percent <- sboptions("percent",FALSE); opt.digits.pval <- sboptions("digits.pval",1)
  ## which pvalue to use
  bp <- "k.1"
  auk <- paste("k.",k,sep="")
  sik <- paste("sk.",k,sep="")
  
  ## table for trees
  a0 <- table.relltest(x$trees)
  a1 <- table.summary.scalebootv(summary(x$trees,k=c(1,k)))$out
  a1[a1[,"hypothesis"] == "null",sik] <- ""
  a <- cbind(a1,a0)
  out.tree <- a[,c("stat","shtest",bp,auk,sik,"beta0","beta1")]
  if(!is.null(x$treename)) {
    out.tree <- cbind(out.tree, x$treename); colnames(out.tree)[ncol(out.tree)] <- "tree"
  }
  out.tree <- cbind(out.tree, sapply(x$tree2edge,function(a) paste("E",sort(a),sep="",collapse=",")))
  colnames(out.tree)[ncol(out.tree)] <- "edge"
  
  ## table for edges
  a1 <- table.summary.scalebootv(summary(x$edges,k=c(1,k)))$out
  a1[a1[,"hypothesis"] == "null",sik] <- ""
  out.edge <- a1[,c(bp,auk,sik,"beta0","beta1")]
  if(!is.null(x$edgename))   {
    out.edge <- cbind(out.edge, x$edgename); colnames(out.edge)[ncol(out.edge)] <- "edge"
  }
  out.edge <- cbind(out.edge, sapply(x$edge2tree,function(a) paste("T",sort(a),sep="",collapse=",")))
  colnames(out.edge)[ncol(out.edge)] <- "tree"
  
  sboptions("percent", opt.percent); sboptions("digits.pval", opt.digits.pval)
  list(tree=out.tree,edge=out.edge,taxa=x$taxaname)
}


