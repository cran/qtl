######################################################################
#
# find.errors.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# Oct, 2001; July, 2001; Apr, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: find.errors, plot.errors
#
######################################################################

######################################################################
#
# find.errors: uses the results of argmax.geno to identify possible
#              genotyping errors
#
######################################################################

find.errors <-
function(cross, chr, error.prob=0.01,
         map.function=c("haldane","kosambi","c-f"),msg=TRUE)
{
  map.function <- match.arg(map.function)
  if(!missing(chr)) cross <- pull.chr(cross,chr)

  # remove chromosomes with < 2 markers
  n.mar <- nmar(cross)
  cross <- pull.chr(cross,names(n.mar)[n.mar >= 2])

  chr <- ind <- mar <- obs <- am <- NULL

  for(i in 1:nchr(cross)) {
    if(is.na(match("argmax",names(cross$geno[[i]])))) { # haven't yet run argmax
      warning("First running argmax.geno.")
      cross <- argmax.geno(cross,error.prob=error.prob,
                           map.function=map.function)
    }
    else { # pull out only marker loci from argmax information
      u <- grep("^loc\-*[0-9]+",colnames(cross$geno[[i]]$argmax))
      if(any(u))
        cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-u]
    }

    X <- cross$geno[[i]]$data
    Y <- cross$geno[[i]]$argmax

    if(class(cross)[1]=="bc") {
      wh <- (!is.na(X) & X != Y)
      geno <- c("AA","AB")
    }
    else if(class(cross)[1]=="f2") {
      wh <- (!is.na(X) & X != Y & !(X == 4 & (Y==1 | Y==2)) & !(X==5 & (Y==2 | Y==3)))
      geno <- c("AA","AB","BB","AA/AB","AB/BB")
    }
    else if(class(cross)[1]=="4way") {
      wh <- (!is.na(X) & X != Y & !(X == 5 & (Y==1 | Y==3)) & !(X==6 & (Y==2 | Y==4)) &
             !(X==7 & (Y==1 | Y==2)) & !(X==8 & (Y==3 | Y==4)) &
             !(X==9 & (Y==1 | Y==4)) & !(X==10 & (Y==2 | Y==3)))
      geno <- c("AC","BC","AD","BD","AC/AD","BC/BD","AC/BC","AD/BD","AC/BD","AD/BC")
    }
    else
      stop(paste("find.errors not available for cross type", class(cross)[1],"."))
      
    if(any(wh)) {
      chr <- c(chr,rep(names(cross$geno)[i],sum(wh)))
      ind <- c(ind,row(X)[wh])
      mar <- c(mar,colnames(X)[col(X)[wh]]) 
      obs <- c(obs,geno[X[wh]])
      am <- c(am,geno[Y[wh]])
    }
                    
  }
  if(is.null(am)) {
    if(msg) cat("\tNo errors found.\n");
    return(invisible(NULL))
  }

  errs <- data.frame(chr=chr,ind=ind,marker=mar,
                     geno=obs,argmax.geno=am)[order(ind),]
  rownames(errs) <- 1:nrow(errs)
  errs
}

  

plot.errors <-
function(x, chr, ...)
{
  cross <- x
  if(!missing(chr)) cross <- pull.chr(cross,chr)

  # remove chromosomes with < 2 markers
  n.mar <- nmar(cross)
  cross <- pull.chr(cross,names(n.mar)[n.mar >= 2])

  errs <- NULL

  for(i in 1:nchr(cross)) {

    if(is.na(match("argmax",names(cross$geno[[i]])))) { # haven't yet run argmax 
      warning("First running argmax.geno.")
      cross <- argmax.geno(cross)
    }
    else { # pull out only marker loci from argmax information
      u <- grep("^loc\-*[0-9]+",colnames(cross$geno[[i]]$argmax))
      if(any(u))
        cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-u]
    }

    X <- cross$geno[[i]]$data
    Y <- cross$geno[[i]]$argmax

    if(class(cross)[1]=="bc") 
      wh <- (!is.na(X) & X != Y)
    else if(class(cross)[1]=="f2") 
      wh <- (!is.na(X) & X != Y & !(X == 4 & (Y==1 | Y==2)) & !(X==5 & (Y==2 | Y==3)))
    else if(class(cross)[1]=="4way")
      wh <- (!is.na(X) & X != Y & !(X == 5 & (Y==1 | Y==3)) & !(X==6 & (Y==2 | Y==4)) &
             !(X==7 & (Y==1 | Y==2)) & !(X==8 & (Y==3 | Y==4)) &
             !(X==9 & (Y==1 | Y==4)) & !(X==10 & (Y==2 | Y==3)))
    else
      stop(paste("find.errors not available for cross type", class(cross)[1],"."))

    nc <- ncol(X); nr <- nrow(X)
    X <- matrix(1,ncol=nc,nrow=nr)
    if(any(wh)) X[wh] <- NA
    cross$geno[[i]]$data <- X
  }

  plot.missing(cross,main="Likely genotyping errors")
}
  
         
# end of find.errors.R
