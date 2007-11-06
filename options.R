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
###
### OPTIONS
###

## model names
sbmodelnames <- function(m=3,one.surface=TRUE,two.surfaces=FALSE,
                         poly,sing,poa,pob,sia,sib,sphe){
  if(missing(poly)) poly <- if(one.surface) m else 0
  if(missing(sing)) sing <- if(one.surface) m else 0
  if(missing(poa)) poa <- if(two.surfaces) m else 0
  if(missing(pob)) pob <- if(two.surfaces) m else 0
  if(missing(sia)) sia <- if(two.surfaces) m else 0
  if(missing(sib)) sib <- if(two.surfaces) m else 0
  if(missing(sphe)) sphe <- 0

  ## one surface models
  polyk <- if(poly>=1) paste("poly",1:poly,sep=".") else NULL
  singk <- if(sing>=3) paste("sing",3:sing,sep=".") else NULL

  ## two surfaces models (polynomial-difference)
  poa1k <- if(poa>=2) paste("poa1",2:poa,sep=".") else NULL
  poa2k <- if(poa>=4) paste("poa2",4:poa,sep=".") else NULL
  pob1k <- if(pob>=2) paste("pob1",2:pob,sep=".") else NULL
  pob2k <- if(pob>=4) paste("pob2",4:pob,sep=".") else NULL

  ## two surfaces models (singular-difference)
  sia1k <- if(sia>=4) paste("sia1",4:sia,sep=".") else NULL
  sia2k <- if(sia>=5) paste("sia2",5:sia,sep=".") else NULL
  sib1k <- if(sib>=4) paste("sib1",4:sib,sep=".") else NULL
  sib2k <- if(sib>=5) paste("sib2",5:sib,sep=".") else NULL

  ## specialized models
  sphek <- if(sphe==3) paste("sphe",sphe,sep=".") else NULL

  ## output model names
  c(polyk,singk,poa1k,poa2k,pob1k,pob2k,sia1k,sia2k,sib1k,sib2k,sphek)
}

## default options (local variable)

.onLoad <- function(libname, pkgname) {
  op <- list(
       ## sbfit
       models = sbmodelnames(), # default models for fitting
       control = list(reltol=1e-14),  # for optim
       method = NULL,  # for optim
       mag.poly = c(1,0.1,0.01,0.001),  # mag factor for par in poly model
       mag.sing = c(1,0.1,0.01,0.001),  # mag factor for par in sing model
       mag1.sing = 0.1,  # mag factor for singularity parameter
       mag1.poa = c(0.1,0.1),  # mag factor for poa models
       mag1.sia = c(0.1,0.1),  # mag factor for sia models
       mag.sphe = c(1.0,1.0,0.01), # mag factor for spherical model
       percent = TRUE, # print p-values in percent
       digits.pval = 2, # significant digits for pvalue
       digits.coef = 4, # significant digits for coefficients
       th.aicw = 1e-4,  # threshold for akaike weights

       ## sbconf; only experimental...
       probs0=c(0.001,0.01,0.1,0.9,0.99,0.999), # initial grid
       prob0=0.5, # for starting value
       tol.mono=0.1,  # for monotonicity checking
       tol.conv=0.01, # for convergence
       tol.relconv=0.5, # tolerance with respect to s.e. 
       max.loop=100, # max iterations

       ## for debug
       debug=FALSE
       )
  options(scaleboot=op)
}


## to show and change options
sboptions <- function(x,value) {
  op <- getOption("scaleboot")
  if(missing(x)) return(op)
  if(length(x)==1) y <- op[[x]]
  else y <- op[x]
  if(!missing(value)) {
    op[[x]] <- value
    options(scaleboot=op)
  }
  y
}
