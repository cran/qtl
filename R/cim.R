######################################################################
#
# cim.R
#
# copyright (c) 2007, Karl W Broman
# 
# last modified Jan, 2007
# first written Jan, 2007
#
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: cim, forwsel
#
######################################################################

######################################################################
# CIM by first doing forward selection at the markers (with filled-in
# data) to a fixed number of markers, followed by interval mapping with
# the selected markers as covariates, dropping marker covariates if
# they are within some fixed window of the location under test.
######################################################################
cim <-
function(cross, pheno.col=1, n.marcovar=3, window=10,
         method=c("em", "imp", "hk", "ehk"),
         imp.method=c("imp", "argmax"), error.prob=0.0001,
         map.function=c("haldane", "kosambi", "c-v", "morgan"),
         n.perm)
{
  method <- match.arg(method)
  imp.method <- match.arg(imp.method)
  map.function <- match.arg(map.function)

  y <- cross$pheno[,pheno.col]
  if(any(is.na(y))) {
    cross <- subset(cross, ind=(!is.na(y)))
    y <- y[!is.na(y)]
  }

  if(!missing(n.perm) && n.perm > 0) {
    results <- rep(NA, n.perm)
    for(i in 1:n.perm) {
      o <- sample(length(y))
      y <- y[o]
      cross$pheno[,pheno.col] <- y
      temp <- cim(cross, pheno.col=pheno.col, n.marcovar=n.marcovar,
                  window=window, method=method, imp.method=imp.method,
                  error.prob=error.prob, map.function=map.function)
      results[i] <- max(temp[,3], na.rm=TRUE)
    }
    class(results) <- "scanoneperm"
    return(results)
  }

  window <- window/2 # window specifies twice the distance between marker and test position

  g <- pull.geno(fill.geno(cross, method=imp.method, error.prob=error.prob,
                           map.function=map.function))

  out.forw <- forwsel(g, y, n.marcovar)
  mar <- colnames(g)[out.forw]
  chrpos <- find.markerpos(cross, mar)

  chrpos

  ac <- g[,mar,drop=FALSE]

  firstscan <- scanone(cross, pheno.col=pheno.col, addcovar=ac,
                       method=method)

  # scan again, dropping one marker covariate at a time
  for(i in seq(along=mar)) {
    temp <- scanone(cross, pheno.col=pheno.col, addcovar=ac[,-i],
                               method=method, chr=chrpos[i,1])
    wh1 <- (firstscan[,1]==chrpos[i,1] & firstscan[,2] >= chrpos[i,2]-window &
            firstscan[,2] <= chrpos[i,2]+window)
    wh2 <- (temp[,2] >= chrpos[i,2]-window & temp[,2] <= chrpos[i,2] + window)

    firstscan[wh1,3] <- temp[wh2,3]
  }

  attr(firstscan, "marker.covar") <- mar
  attr(firstscan, "marker.covar.pos") <- chrpos
  
  u <- table(chrpos[,1])
  if(any(u>1)) {
    u <- names(u)[u>1]
    wh <- which(chrpos[,1]==u)
    pos <- chrpos[wh,2]
    d <- diff(pos)
    if(any(d <= window)) {
      w <- which(d <= window)
      for(i in seq(along=w)) {
        tempac <- ac[,-wh[c(w[i],w[i]+1)]]
        thepos <- pos[c(w[i],w[i]+1)]
        temp <- scanone(cross, pheno.col=pheno.col, addcovar=tempac,
                        method=method, chr=u)
        wh1 <- (firstscan[,1]==u & firstscan[,2] >= thepos[2]-window &
                firstscan[,2] <= thepos[1]+window)
        wh2 <- (temp[,2] >= thepos[2]-window & temp[,2] <= thepos[1]-window)
        firstscan[wh1,3] <- temp[wh1,3]
      }
    }
        
  }
  firstscan
}

######################################################################
# Simple forward selection to a fixed number of covariates
#
# x = matrix of covariates
# y = outcome
# maxsize = maximum size of model
#
# output: indices of chosen covariates [1, 2, ..., ncol(x)]
######################################################################
forwsel <-
function(x, y, maxsize=7)
{
  if(length(y) != nrow(x))
    stop("Need length(y) == nrow(x).")

  if(maxsize < 0 || maxsize > ncol(x))
    stop("Need maxsize between 1 and ncol(x).")

  out <- .C("R_forwsel",
            as.integer(nrow(x)),
            as.integer(ncol(x)),
            as.double(x),
            as.double(y),
            as.integer(maxsize),
            chosen=as.integer(rep(0,maxsize)),
            rss=as.double(rep(0,maxsize)),
            PACKAGE="qtl")

  out$chosen+1
}

# end of cim.R
