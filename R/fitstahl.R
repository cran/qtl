######################################################################
#
# fitstahl.R
#
# copyright (c) 2006, Karl W Broman, Johns Hopkins University
# last modified Aug, 2006
# first written Aug, 2006
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: fitstahl, fitstahl.estp, fitstahl.este, fitstahl.estpe,
#           fitstahl.estp.sub, fitstahl.este.sub, fitstahl.estpe.sub
#
######################################################################

######################################################################
# fitstahl
#
# Fit the Stahl model for crossover interference (or the chi-square
# model, which is a special case)
#
# cross: the cross object
#
# chr:   Chromosome(s) to use; if unspecified, pooled estimates for
#        all chromosomes are obtained
#
# m:     Interference parameter (a non-negative integer); if unspecified,
#        this is estimated
#
# p:     The proportion of chiasmata coming from the no interference
#        mechanism in the Stahl model (0 <= p <= 1).  p=0 gives the
#        chi-square model.  If unspecified, this is estimated
#
# error.prob  The genotyping error probability.  If = NULL, it is
#             estimated
#
# NOTE: If m, p, error.prob are all specified, they can be vectors
#       or have length 1; any with length > 1 should all have the
#       same length
#
# maxit  Maximum number of iterations
#
# tol    Tolerance for convergence
#
# maxm   Maximum value of m to consider, if m is unspecified.
#
######################################################################

fitstahl <-
function(cross, chr, m, p, error.prob=0.0001, maxit=4000, tol=1e-4,
         maxm=15, verbose=TRUE)
{
  if(length(class(cross)) < 2 || class(cross)[2] != "cross")
    stop("Input should have class \"cross\".")

  if(class(cross)[1] != "bc")
    stop("fitstahl only working for backcrosses.")

  if(!missing(chr)) cross <- subset(cross, chr)

  if(!missing(m)) { # m was specified
    if(!missing(p)) { # p was specified

      # m, p, error.prob all specified
      if(!is.null(error.prob)) { 

        n <- c(length(m), length(p), length(error.prob))
        mn <- max(n)
        if(mn > 1) {
          if(any(n > 1 & n < mn))
            stop("Any m, p, error.prob with length > 1 must have same length")
          if(length(m) == 1) m <- rep(m, mn)
          if(length(p) == 1) p <- rep(p, mn)
          if(length(error.prob) == 1) error.prob <- rep(error.prob, mn)
        }

        result <- matrix(ncol=4, nrow=mn)
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        for(i in 1:mn) {
          result[i,] <- c(m[i], p[i], error.prob[i],
                          fitstahl.estp.sub(p[i], cross, error.prob[i],
                                            m[i], tol, maxit))
          if(verbose) cat(i,result[i,], "\n")
        }
      }

      # m,p specified but error.prob wasn't
      else {
        if(length(m) == 1)
          m <- rep(m, length(p))
        else {
          if(length(p) == 1)
            p <- rep(p, length(m))
          else if(length(m) != length(p))
            stop("Any m, p, error.prob with length > 1 must have same length")
        }
        mn <- length(m)

        result <- matrix(ncol=4, nrow=mn)
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        for(i in 1:mn) {
          temp <- fitstahl.este(cross, m[i], p[i], tol, maxit)
          result[i,] <- c(m[i], p[i], temp$est, temp$loglik)
          if(verbose) cat(i, result[i,], "\n")
        }
      }
    }
    else {
      # m, error.prob specified; p unspecified
      if(!is.null(error.prob)) { 
        if(length(m) == 1)
          m <- rep(m, length(error.prob))
        else {
          if(length(error.prob) == 1)
            error.prob <- rep(error.prob, length(m))
          else if(length(m) != length(error.prob))
            stop("Any m, p, error.prob with length > 1 must have same length")
        }
        mn <- length(m)

        result <- matrix(ncol=4, nrow=mn)
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        for(i in 1:mn) {
          temp <- fitstahl.estp(cross, error.prob[i], m[i], tol, maxit)
          result[i,] <- c(m[i], temp$est, error.prob[i], temp$loglik)
          if(verbose) cat(i, result[i,], "\n")
        }
      }

      # only m specified
      else { 
        result <- matrix(ncol=4, nrow=length(m))
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        for(i in 1:length(m)) {
          temp <- fitstahl.estpe(cross, m[i], tol, maxit)
          result[i,] <- c(m[i], temp$est, temp$loglik)
          if(verbose) cat(i, result[i,], "\n")
        }
      }

    }
  }
  else { # m unspecified
    if(!missing(p)) { # p was specified

      # p, e specified
      if(!is.null(error.prob)) { 
        if(length(p) == 1)
          p <- rep(p, length(error.prob))
        else {
          if(length(error.prob) == 1)
            error.prob <- rep(error.prob, length(p))
          else if(length(p) != length(error.prob))
            stop("Any m, p, error.prob with length > 1 must have same length")
        }
        mn <- length(p)

        result <- matrix(ncol=4, nrow=mn)
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        for(i in 1:mn) {
          # fit the case m=0
          maxll <- fitstahl.estp.sub(p[i], cross, error.prob[i], 0, tol, maxit)
          themax <- 0
          if(verbose) cat(i, 0, maxll, "\n")
          for(j in 1:maxm) {
            curll <- fitstahl.estp.sub(p[i], cross, error.prob[i], j, tol, maxit)
            if(verbose) cat(i, j, curll, "\n")
            if(curll < maxll) break
            if(curll > maxll) {
              maxll <- curll
              themax <- j
            }
          }
          result[i,] <- c(themax, p[i], error.prob[i], maxll)
          if(verbose) cat(i, result[i,], "\n")
        }
      }

      # only p specified
      else { 
        mn <- length(p)

        result <- matrix(ncol=4, nrow=mn)
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        for(i in 1:mn) {
          # fit the case m=0
          temp <- fitstahl.este(cross, 0, p[i], tol, maxit)
          maxll <- temp$loglik
          themax <- 0
          themaxe <- temp$est
          if(verbose) cat(i, 0, maxll, "\n")
          for(j in 1:maxm) {
            temp <- fitstahl.este(cross, j, p[i], tol, maxit)
            curll <- temp$loglik
            if(verbose) cat(i, j, curll, "\n")
            if(curll < maxll) break
            if(curll > maxll) {
              maxll <- curll
              themax <- j
              themaxe <- temp$est
            }
          }
          result[i,] <- c(themax, p[i], themaxe, maxll)
          if(verbose) cat(i, result[i,], "\n")
        }
      }
    }
    else {
      if(!is.null(error.prob)) { # error.prob specified; p unspecified
        mn <- length(error.prob)

        result <- matrix(ncol=4, nrow=mn)
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        for(i in 1:mn) {
          # fit the case m=0 (in which case p doesn't matter)
          maxll <- fitstahl.estp.sub(0, cross, error.prob[i], 0, tol, maxit)
          themax <- 0
          themaxp <- 0
          if(verbose) cat(i, 0, maxll, "\n")
          for(j in 1:maxm) {
            temp <- fitstahl.estp(cross, error.prob[i], j, tol, maxit)
            curll <- temp$loglik
            if(verbose) cat(i, j, curll, "\n")
            if(curll < maxll) break
            if(curll > maxll) {
              maxll <- curll
              themax <- j
              themaxp <- temp$est
            }
          }
          result[i,] <- c(themax, themaxp, error.prob[i], maxll)
          if(verbose) cat(i, result[i,], "\n")
        }
      }

      else { # nothing specified

        result <- matrix(ncol=4, nrow=1)
        colnames(result) <- c("m", "p", "error.prob", "loglik")
        # fit the case m=0 (in which case p doesn't matter)
        temp <- fitstahl.este(cross, 0, 0, tol, maxit)
        maxll <- temp$loglik
        themax <- 0
        themaxest <- c(0, temp$est) 
        if(verbose) cat(i, 0, maxll, "\n")
        for(j in 1:maxm) {
          temp <- fitstahl.estpe(cross, j, tol, maxit)
          curll <- temp$loglik
          if(verbose) cat(i, j, curll, "\n")
          if(curll < maxll) break
          if(curll > maxll) {
            maxll <- curll
            themax <- j
            themaxest <- temp$est
          }
        }
        result[1,] <- c(themax, themaxest, maxll)
        if(verbose) cat(1, result[1,], "\n")
      }
    }
  }
  as.data.frame(result)
}
    


######################################################################
# fitstahl.estp:  estimate p for fixed m and error.prob
######################################################################
fitstahl.estp <-
function(cross, error.prob=0.0001, m=0, tol=1e-4, maxit=4000)
{
  out <- optimize(fitstahl.estp.sub, interval=c(0,1), maximum=TRUE,
                  cross=cross, m=m, error.prob=error.prob, thetol=tol,
                  maxit=maxit)

  # make sure we fit p=0
  temp <- fitstahl.estp.sub(0, cross, error.prob, m, tol, maxit)
  if(temp >= out[[2]]) {
    est <- 0
    loglik <- temp
  }
  else {
    est <- out[[1]]
    loglik <- out[[2]]
  }
  list(est=est, loglik=loglik)
}

fitstahl.estp.sub <-
function(p, cross, error.prob=0.0001, m=0, thetol=1e-4, maxit=4000)
  sum(sapply(est.map(cross, error.prob=error.prob, m=m, p=p, tol=thetol,
                     maxit=maxit), function(a) attr(a, "loglik")))

######################################################################
# fitstahl.este: estimate error.prob for fixed m and p
######################################################################
fitstahl.este <-
function(cross, m=0, p=0, tol=1e-4, maxit=4000)
{
  out <- optimize(fitstahl.este.sub, interval=c(0,1), maximum=TRUE,
                  cross=cross, m=m, p=p, thetol=tol, maxit=maxit)

  # make sure we fit error.prob=0
  temp <- fitstahl.este.sub(0, cross, m, p, tol, maxit)
  if(temp >= out[[2]]) {
    est <- 0
    loglik <- temp
  }
  else {
    est <- out[[1]]
    loglik <- out[[2]]
  }

  list(est=est, loglik=loglik)
}

fitstahl.este.sub <-
function(error.prob, cross, m=0, p=0, thetol=1e-4, maxit=4000)
  sum(sapply(est.map(cross, error.prob=error.prob, m=m, p=p, tol=thetol,
                     maxit=maxit), function(a) attr(a, "loglik")))


######################################################################
# fitstahl.estpe: estimate p and error.prob for fixed m 
######################################################################
fitstahl.estpe <-
function(cross, m=0, tol=1e-4, maxit=4000)
{
  out <- optim(c(0.1,0.01), fitstahl.estpe.sub, method="L-BFGS-B",
               lower=c(0,0), upper=c(1,1),
               control=list(fnscale=-1, maxit=maxit, reltol=tol),
               cross=cross, m=m, thetol=tol, maxit=maxit)

  if(out$convergence !=0) 
    warning(" Didn't converge.")
  print(out$counts)

  # make sure we fit the case p=0
  temp <- fitstahl.estp(cross, m, 0, tol, maxit)
  templl <- attr(temp, "loglik")
  if(templl >= out$value) {
    est <- c(0, temp) 
    loglik <- templl
  }
  else {
    est <- out$par
    loglik <- out$value
  }
  names(est) <- c("p", "error.prob")
  list(est=est, loglik=loglik)
}

fitstahl.estpe.sub <-
function(x, cross, m=0, thetol=1e-4, maxit=4000)
  sum(sapply(est.map(cross, error.prob=x[2], m=m, p=x[1], tol=thetol,
                     maxit=maxit), function(a) attr(a, "loglik")))


# end of fitstahl.R
