######################################################################
#
# est.map.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# last modified Sept, 2001
# first written Apr, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: est.map
#
######################################################################

######################################################################
#
# est.map: re-estimate the genetic map for an experimental cross
#
######################################################################

est.map <- 
function(cross, error.prob=0, map.function=c("haldane","kosambi","c-f"),
         maxit=4000, tol=1e-4, sex.sp=TRUE, trace=FALSE)
{

  # map function
  map.function <- match.arg(map.function)
  if(map.function=="kosambi") {
    mf <- mf.k; imf <- imf.k
  }
  else if(map.function=="c-f") {
    mf <- mf.cf; imf <- imf.cf
  }
  else {
    mf <- mf.h; imf <- imf.h
  }

  # don't let error.prob be exactly zero, just in case
  if(error.prob < 1e-50) error.prob <- 1e-50

  n.ind <- nind(cross)
  n.mar <- nmar(cross)
  n.chr <- nchr(cross)

  newmap <- vector("list",n.chr)
  names(newmap) <- names(cross$geno)

  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {

    if(n.mar[i] < 2) {
      newmap[[i]] <- cross$geno[[i]]$map
      next
    }

    # which type of cross is this?
    if(class(cross)[1] == "f2") {
      one.map <- TRUE
      if(class(cross$geno[[i]]) == "A") # autosomal
        cfunc <- "est_map_f2"
      else                              # X chromsome 
        cfunc <- "est_map_bc"
    }
    else if(class(cross)[1] == "bc") {
      one.map <- TRUE
      cfunc <- "est_map_bc"
    }
    else if(class(cross)[1] == "4way") {
      one.map <- FALSE
      cfunc <- "est_map_4way"
    }
    else if(class(cross)[1] == "f2ss") {
      one.map <- FALSE
      cfunc <- "est_map_f2ss"
    }
    else stop(paste("est.map not available for cross type",
                    class(cross)[1], "."))

    # genotype data
    gen <- cross$geno[[i]]$data
    gen[is.na(gen)] <- 0
    
    # recombination fractions
    if(one.map) {
      # recombination fractions
      rf <- mf(diff(cross$geno[[i]]$map))
      rf[rf < 1e-14] <- 1e-14
    }
    else {
      rf <- mf(diff(cross$geno[[i]]$map[1,]))
      rf[rf < 1e-14] <- 1e-14
      rf2 <- mf(diff(cross$geno[[i]]$map[2,]))
      rf2[rf2 < 1e-14] <- 1e-14
      if(!sex.sp && class(cross$geno[[i]])=="X")
        temp.sex.sp <- TRUE
      else temp.sex.sp <- sex.sp
    }


    if(trace) cat(paste("Chr ", names(cross$geno)[i], ":\n",sep="")) 

    # call the C function
    if(one.map) {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.mar[i]),      # number of markers
              as.integer(gen),           # genotype data
              rf=as.double(rf),          # recombination fractions
              as.double(error.prob),     
              loglik=as.double(0),       # log likelihood
              as.integer(maxit),
              as.double(tol),
              as.integer(trace),
              PACKAGE="qtl")

      newmap[[i]] <- cumsum(c(min(cross$geno[[i]]$map),imf(z$rf)))
      names(newmap[[i]]) <- names(cross$geno[[i]]$map)
      attr(newmap[[i]],"loglik") <- z$loglik
    }
    else {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.mar[i]),      # number of markers
              as.integer(gen),           # genotype data
              rf=as.double(rf),          # recombination fractions
              rf2=as.double(rf2),        # recombination fractions
              as.double(error.prob),
              loglik=as.double(0),       # log likelihood
              as.integer(maxit),
              as.double(tol),
              as.integer(temp.sex.sp),
              as.integer(trace),
              PACKAGE="qtl")
              
      if(!temp.sex.sp) z$rf2 <- z$rf

      newmap[[i]] <- rbind(cumsum(c(min(cross$geno[[i]]$map[1,]),imf(z$rf))),
                           cumsum(c(min(cross$geno[[i]]$map[2,]),imf(z$rf2))))
      dimnames(newmap[[i]]) <- dimnames(cross$geno[[i]]$map)
      attr(newmap[[i]],"loglik") <- z$loglik
    }

  } # end loop over chromosomes

  class(newmap) <- "map"
  newmap
}

# end of est.map.R
