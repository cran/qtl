######################################################################
#
# ripple.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# Oct, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: ripple, summary.ripple
#           ripple.perm1, ripple.perm2, ripple.perm.sub
#
######################################################################

######################################################################
#
# ripple: Check marker orders for a given chromosome, comparing all
#         possible permutations of a sliding window of markers
#
######################################################################

ripple <-
function(cross, chr, window=4, error.prob=0,
         map.function=c("haldane","kosambi","c-f"),
         maxit=1000,tol=1e-5,sex.sp=TRUE)
{
  # pull out relevant chromosome
  if(length(chr) > 1)
    stop("ripple only works for one chromosome at a time.")
  cross <- pull.chr(cross,chr)
  chr.name <- names(cross$geno)[1]

  map.function <- match.arg(map.function)

  # get marker orders to test
  n.mar <- totmar(cross)
  if(n.mar <= window) # look at all possible orders
    orders <- ripple.perm2(n.mar)
  else { 
    temp <- ripple.perm1(window)
    n <- nrow(temp)
    orders <- cbind(temp,matrix(rep((window+1):n.mar,n),
                                  byrow=TRUE,ncol=n.mar-window))
    for(i in 2:(n.mar-window+1)) {
      left <- matrix(rep(1:(i-1),n),byrow=TRUE,ncol=i-1)
      if(i < n.mar-window+1)
        right <- matrix(rep((i+window):n.mar,n),byrow=TRUE,ncol=n.mar-window-i+1)
      else
        right <- NULL
      orders <- rbind(orders,cbind(left,temp+i-1,right))
    }
    # keep only distinct orders
    orders <- as.numeric(unlist(strsplit(unique(apply(orders,1,paste,collapse=":")),":")))
    orders <- matrix(orders,ncol=n.mar,byrow=TRUE)
  }

  m <- seq(0,by=5,length=n.mar)
  temcross <- cross
  if(is.matrix(cross$geno[[1]]$map)) 
    temcross$geno[[1]]$map <- rbind(m,m)
  else temcross$geno[[1]]$map <- m

  # calculate log likelihoods (and est'd chr length) for each marker order
  n.orders <- nrow(orders)
  loglik <- 1:n.orders
  chrlen <- 1:n.orders

  # how often to print information about current order being considered
  if(n.orders > 49) print.by <- 10
  else if(n.orders > 14) print.by <- 5
  else print.by <- 2

  for(i in 1:n.orders) {
    if(i==1) cat("  ", n.orders,"total orders\n")
    if((i %/% print.by)*print.by == i) cat("    --Order", i, "\n")
    temcross$geno[[1]]$data <- cross$geno[[1]]$data[,orders[i,]]
    newmap <- est.map(temcross,error.prob,map.function,maxit,tol,sex.sp)
    loglik[i] <- attr(newmap[[1]],"loglik")
    chrlen[i] <- newmap[[1]][n.mar]
  }

  # re-scale log likelihoods and convert to lods
  loglik <- (loglik - loglik[1])/log(10)
  
  # sort orders by lod
  o <- rev(order(loglik[-1])+1)

  orders <- cbind(orders,lod=loglik,chrlen)[c(1,o),]
  class(orders) <- c("ripple","matrix")
  attr(orders,"chr") <- chr.name
  attr(orders,"window") <- window
  attr(orders,"error.prob") <- error.prob
  orders
}

######################################################################
#
# summary.ripple: print top results from ripple().  We do this so
#                 that we can return *all* results but allow easy
#                 view of only the important ones
#
######################################################################

summary.ripple <-
function(object,lod.cutoff=2,...)
{
  n <- ncol(object)
  if(!any(object[-1,n-1] > -lod.cutoff)) object <- object[1:2,]
  else # make sure first row is included
    object <- rbind(object[1,],object[-1,][object[-1,n-1] > -lod.cutoff,])

  class(object) <- c("summary.ripple","matrix")
  object
}

######################################################################
#
# ripple.perm1: Utility function for ripple().  Returns all possible
#               permutations of {1, 2, ..., n}
#
######################################################################

ripple.perm1 <-  
function(n)
{
  if(n == 1) return(rbind(1))
  o <- rbind(c(n-1,n),c(n,n-1))
  if(n > 2)
    for(i in (n-2):1)
      o <- ripple.perm.sub(i,o)
  dimnames(o) <- NULL
  o
}

######################################################################
#
# ripple.perm2: Utility function for ripple().  Returns all possible
#               permutations of {1, 2, ..., n}, up to orientation of
#               the entire group
#
######################################################################

ripple.perm2 <- 
function(n)
{
  if(n < 3) return(rbind(1:n))
  o <- rbind(c(n-2,n-1,n),c(n-1,n-2,n),c(n-1,n,n-2))
  if(n > 3)
    for(i in (n-3):1)
      o <- ripple.perm.sub(i,o)
  dimnames(o) <- NULL
  o
}

######################################################################
#
# ripple.perm.sub: Subroutine used for ripple().  I'm too tired to
#                  explain.
#
######################################################################

ripple.perm.sub <-
function(x,mat)
{
  res <- cbind(x,mat)
  if(ncol(mat) > 1) {
    for(i in 1:ncol(mat))
        res <- rbind(res,cbind(mat[,1:i],x,mat[,-(1:i)]))
  }
  res
}

# end of ripple.R
