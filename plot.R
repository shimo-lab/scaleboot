##
##  scaleboot: R package for multiscale bootstrap
##  Copyright (C) 2006-2007 Hidetoshi Shimodaira
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
### MAIN: PLOT

##
## mb : output of sbfit
## xval : x-axis variable
## yval : y-axis variable
## xlog : plot in log for x-axis
## ylog : plot in log for y-axis
##

plot.scaleboot <- function(x,
                       models=names(x$fi),
                       k=NULL,s=NULL,sp=NULL,lambda=NULL,
                       ## parameters for x-y axises
                       xval=c("square","inverse","sigma"),
                       yval=c("psi","zvalue","pvalue"),
                       xlab=NULL,ylab=NULL,
                       log.xy="",
                       xlim=NULL,ylim=NULL,
                       add=F,
                       length.x=300,
                       main=NULL,
                       ## lines
                       col=1:6,lty=1:5,lwd=par("lwd"),
                       ## extrapolation
                       ex.pch=2:7,
                       ## points
                       pch=1,cex=1,pt.col=col[1],pt.lwd=lwd[1],
                       ## legend
                       legend.x=NULL, inset=0.1,
                       ...
                       ) {

  ## check log
  if(length(log.xy) && log.xy != "") a <- strsplit(log.xy,"")[[1]] else a <- ""
  xlog <- "x" %in% a
  ylog <- "y" %in% a

  ## specify models
  if(is.numeric(models)) models <- names(x$fi)[models]
  else if(length(models)==1 && models=="best") models <- x$best$model
  else if(length(models)==1 && models=="average") models <-
                            names(sort(x$average$w,decreasing=TRUE))
  ## check if models are valid
  i <- match(models,names(x$fi))
  models <- names(x$fi)[i[!is.na(i)]]
  
  ## x-axis
  ## s : sigma^2
  ## xfun(s) : x-axis function
  ## xinv(x) : s (inverse of xfun)
  xval <- match.arg(xval)
  if(xval=="sigma") {
    xlab0 <- expression(sigma)
    xfun <- function(s) s^0.5
    xinv <- function(x) x^2
  } else if(xval=="inverse") {
    xlab0 <- expression(1/sigma)
    xfun <- function(s) s^(-0.5)
    xinv <- function(x) x^(-2)
  } else if(xval=="square") {
    xlab0 <- expression(sigma^2)
    xfun <- function(s) s
    xinv <- function(x) x
  } else stop("xval is ",xval)
  if(is.null(xlab)) xlab <- xlab0

  ## y-axis
  ## p : psi-value
  ## a : probability-value
  ## yfun(p,s) : conversion p -> y
  ## yfun2(a,s) : conversion a -> y
  yval <- match.arg(yval)
  if(yval=="zvalue") {
    ylab0 <- expression(z(sigma^2))
    yfun <- function(p,s) p/sqrt(s)
    yfun2 <- function(a,s) -qnorm(a)
  } else if(yval=="pvalue") {
    ylab0 <- expression(alpha(sigma^2))
    yfun <- function(p,s) pnorm(-p/sqrt(s))
    yfun2 <- function(a,s) a
  } else if(yval=="psi") {
    ylab0 <- expression(psi(sigma^2))
    yfun <- function(p,s) p
    yfun2 <- function(a,s) -sqrt(s)*qnorm(a)
  } else stop("yval is ",yval)
  if(is.null(ylab)) ylab <- ylab0

  ## observed points
  sa <- x$sa # sigma squared
  bp <- x$bp # bootstrap probabilities
  ss <- sqrt(sa) # scales
  bs <- -qnorm(bp)*ss # observed psi-value
  by <- yfun(bs,sa) # observed y-axis
  bx <- xfun(sa) # observed x-axis
  u <- is.finite(if(ylog) log(by) else by)
  by.u <- by[u]; bx.u <- bx[u]

  ## extrapolation
  if(!is.null(k) && length(models)>0) {
    fis <- x$fi[models]
    zex <- t(sapply(fis,function(fi) {
      psi <- get(fi$psi)
      sapply(k,function(k0) psi(fi$mag*fi$par,s,k=k0,sp=sp,lambda=lambda))}))
    if(length(models)>1) {
      w <- x$average$w[models]
      w <- w/sum(w)
      pex <- apply(w*pnorm(-zex),2,sum) # p-value
      py <- yfun2(pex,1)
    } else {
      py <- yfun(drop(zex),1)
    }
    px <- rep(xfun(sp),length(py))
#    if(is.null(main)) main <- paste("extrapolation (",
#                                    paste(models,collapse="+"),")",sep="")
    if(is.null(main)) main <- c(paste("extrapolation k=",
                                      paste(k,collapse=","),sep=""),
                                paste(models,collapse="+"))
  } else {
    px <- NULL; py <- NULL
    if(is.null(main)) main <- c("model fitting",
                                paste(models,collapse=","))
  }

  ## xlim and ylim
  if(is.null(xlim)) {
    xlim <- if(length(bx.u)>0) range(bx.u,px) else range(bx,px)
  }
  if(is.null(ylim)) {
    ylim <- if(length(bx.u)>0) range(by.u,py) else c(-1,1)
  }

  ## plot
  if(add) {
    points(bx.u,by.u,pch=pch,cex=cex,col=pt.col,lwd=pt.lwd,...)
  }
  else {
    plot(bx.u,by.u,
         xlab=xlab,ylab=ylab,log=log.xy,xlim=xlim,ylim=ylim,
         pch=pch,cex=cex,col=pt.col,lwd=pt.lwd,main=main,...)
  }
  if(!is.null(k))
    points(c(NA,px),c(NA,py),pch=ex.pch,cex=cex,col=col,lwd=pt.lwd,...)


  ## return value (passed to lines)
  z1 <- list(sa=sa,bp=bp,bx=bx,by=by,use=u,
             px=px,py=py,
             xlim=xlim,ylim=ylim,xlog=xlog,ylog=ylog,
             xfun=xfun,xinv=xinv,yfun=yfun,yfun2=yfun2,xlab=xlab,ylab=ylab,
             length.x=length.x,col=col,lty=lty,lwd=lwd)

  ## lines
  z2 <- lines(x,z1,models=models,k=k,s=s,sp=sp,lambda=lambda)


  ## legend
  if(!is.null(legend.x)) sblegend(legend.x,z=c(z1,z2),inset=inset)

  invisible(c(z1,z2))
}

plot.summary.scaleboot <-
  function(x,models="average",
           k=x$parex$k,s=x$parex$s,sp=x$parex$sp,lambda=x$parex$lambda,
           ...)
  plot.scaleboot(x,models=models,k=k,s=s,sp=sp,lambda=lambda,...)


##
## mb : output of sbfit
## z : output of plot.sbfit
##
lines.scaleboot <- function(x,z,
                        models=names(x$fi),
                        k=NULL,s=NULL,sp=NULL,lambda=NULL,
                        length.x=z$length.x,
                        col=z$col,lty=z$lty,lwd=z$lwd,...
                        ) {
  if(is.null(models)||length(models)==0) return(invisible(NULL))
  
  ## sequence in x-axis
  rx <- z$xlim
  if(z$xlog) rx <- log(rx)
  a <- (rx[2]-rx[1])*0.05
  if(z$xlog || !is.null(k)) rx[1] <- rx[1]-a
  else rx[1] <- if(rx[1]>a*3) rx[1] - a else rx[1]*(1-1/3)
  rx[2] <- rx[2]+a
  xx <- seq(from=rx[1],to=rx[2],length=length.x)
  if(z$xlog) xx <- exp(xx)
  sa <- z$xinv(xx)
  u1 <- sa >= 0  # for psi function
  u2 <- sa <= s  # for extraplolation
  
  if(is.null(k)) {### fitting models (without extrapolations)
    yy <- matrix(0,length(xx),length(models))
    yy[] <- NA
    for(i in seq(along=models)) {
      f <- x$fi[[models[i]]]
      if(!is.null(f)) {
        beta <- f$par*f$mag
        psi <- get(f$psi)
        py <- sapply(sa[u1],function(s) psi(beta,s))
        yy[u1,i] <- z$yfun(py,sa)
      }
    }
    labels <- models
    
  } else { ### fitting models and extrapolation
    ## first, prepare the containers
    yy <- matrix(NA,length(xx),length(k)+1) # for matlines
    yy1 <- matrix(NA,sum(u1),length(models)) # for psi
    yy2 <- array(NA,dim=c(sum(u2),length(k),length(models))) # for ext
    ## compute psi values for models
    for(i in seq(along=models)) {
      f <- x$fi[[models[i]]]
      beta <- f$par*f$mag
      psi <- get(f$psi)
      yy1[,i] <- sapply(sa[u1],function(s0) psi(beta,s0))
      for(j in seq(along=k)) {
        yy2[,j,i] <- sapply(sa[u2],function(s1)
                            psi(beta,s=s,k=k[j],sp=s1,lambda=lambda))
      }
    }
    if(length(models)>1) { ## then, average models...
      ## akaike weights
      w <- x$average$w[models]
      w <- w/sum(w)
      ## averaged fitting model
      yy1 <- pnorm(-yy1/sqrt(sa[u1])) # convert psi -> pval
      y1 <- apply(yy1,1,function(y) sum(w*y)) # averaging
      yy[u1,1] <- z$yfun2(y1,sa[u1])
      ## averaged extrapolation
      yy2 <- pnorm(-yy2) # convert psi -> pval
      y2 <- apply(yy2,1:2,function(y) sum(w*y)) # averaging
      yy[u2,-1] <- z$yfun2(y2,1)
    } else { ## if only one model...
      yy[u1,1] <- z$yfun(yy1,sa[u1]) # fitting
      yy[u2,-1] <- z$yfun(yy2,1) # extrapolation
    }
    labels <- c(models,paste("k",k,sep="."))
  }

  if(!all(is.na(yy))) matlines(xx,yy,col=col,lty=lty,lwd=lwd)
  invisible(list(col=col,lty=lty,lwd=lwd,labels=labels,do.lines=T))
}

##
## x,y : legend position
## z : output of plot.sbfit or plot.sbconf
##

sblegend <- function(x="topright",y=NULL,z,inset=0.1,...) {
  if(length(z)==1) z <- z[[1]]
  if(is.null(z$labels)) return(invisible(NULL))
  lty <- rep(z$lty,length=length(z$labels))
  pch <- rep(z$pch,length=length(z$labels))
  if(is.null(do.lines <- z$do.lines)) do.lines <- F 
  if(is.null(do.points <- z$do.points)) do.points <- F
  if(do.lines && ! do.points) 
    legend(x,y,z$labels,col=z$col,lwd=z$lwd,lty=lty,inset=inset,...)
  else if(!do.lines &&  do.points) 
    legend(x,y,z$labels,col=z$col,pch=pch,inset=inset,...)
  else if(do.lines &&  do.points) 
    legend(x,y,z$labels,col=z$col,lwd=z$lwd,lty=lty,pch=pch,inset=inset,...)
  else invisible(NULL)
}

##
## scalebootv
##

plot.scalebootv <- function(x,models=attr(x,"models"),...) {
  ## preliminary
  n <- length(x)
  m <- n2mfrow(n)
  s <- matrix(c(1:n,rep(0,m[1]*m[2]-n)),m[1],m[2],byrow=T)
  def.par <- par(no.readonly = TRUE) # save default
  on.exit(par(def.par))
  layout(s)

  # plots
  z <- NULL
  for(i in seq(length=n)) {
    z[[i]] <- plot(x[[i]],main=names(x)[[i]],models=models,...)
  }
  invisible(z)
}

plot.summary.scalebootv <- function(x, models="average",...)
  plot.scalebootv(x,models=models,...)

