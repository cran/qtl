######################################################################
#
# argmax.geno.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# Sept, 2001; July, 2001; May, 2001; Apr, 2001; Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: argmax.geno
#
######################################################################

######################################################################
#
# argmax.geno: Use Viterbi algorithm to find most likely sequence of
#              underlying genotypes, given observed marker data
#
######################################################################

argmax.geno <-
function(cross, step=0, off.end=0, error.prob=0,
         map.function=c("haldane","kosambi","c-f"))
{

  # map function
  map.function <- match.arg(map.function)
  if(map.function=="kosambi") mf <- mf.k
  else if(map.function=="c-f") mf <- mf.cf
  else mf <- mf.h

  # don't let error.prob be exactly zero, just in case
  if(error.prob < 1e-14) error.prob <- 1e-14

  n.ind <- nind(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)

  # loop over chromosomes
  for(i in 1:n.chr) {

    # which type of cross is this?
    if(class(cross)[1] == "f2") {
      one.map <- TRUE
      if(class(cross$geno[[i]]) == "A") { # autosomal
        cfunc <- "argmax_geno_f2"
      }
      else {                             # X chromsome 
        cfunc <- "argmax_geno_bc"
      }
    }
    else if(class(cross)[1] == "bc") {
      cfunc <- "argmax_geno_bc"
      one.map <- TRUE
    }
    else if(class(cross)[1] == "4way") {
      cfunc <- "argmax_geno_4way"
      one.map <- FALSE
    }
    else if(class(cross)[1] == "f2ss") {
      cfunc <- "argmax_geno_f2ss"
      one.map <- FALSE
    }
    else stop(paste("argmax.geno not available for cross type",
                    class(cross)[1], "."))

    # genotype data
    gen <- cross$geno[[i]]$data
    gen[is.na(gen)] <- 0
    
    # recombination fractions
    if(one.map) {
      # recombination fractions
      map <- create.map(cross$geno[[i]]$map,step,off.end)
      rf <- mf(diff(map))
      rf[rf < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=length(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,names(map))
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
    }
    else {
      map <- create.map(cross$geno[[i]]$map,step,off.end)
      rf <- mf(diff(map[1,]))
      rf[rf < 1e-14] <- 1e-14
      rf2 <- mf(diff(map[2,]))
      rf2[rf2 < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,dimnames(map)[[2]])
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
    }
    if(any(is.na(rf))) { # this occurs when there is only one marker
      rf <- rf2 <- 0
      warning(paste("Only one marker on chr ", names(cross$geno)[i],
                    ": argmax results tenuous.", sep=""))
    }


    # call the C function
    if(one.map) {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(newgen),        # genotype data
              as.double(rf),             # recombination fractions
              as.double(error.prob),     
              argmax=as.integer(newgen), # the output
              PACKAGE="qtl")

      cross$geno[[i]]$argmax <- matrix(z$argmax,ncol=n.pos)
      dimnames(cross$geno[[i]]$argmax) <- list(NULL, names(map))
    }
    else {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(newgen),        # genotype data
              as.double(rf),             # recombination fractions
              as.double(rf2),            # recombination fractions
              as.double(error.prob),      
              argmax=as.integer(newgen), # the output
              PACKAGE="qtl")

      cross$geno[[i]]$argmax <- matrix(z$argmax,ncol=n.pos)
      dimnames(cross$geno[[i]]$argmax) <- list(NULL, colnames(map))
    }

    # attribute set to the error.prob value used, for later
    #     reference
    attr(cross$geno[[i]]$argmax,"error.prob") <- error.prob
    attr(cross$geno[[i]]$argmax,"step") <- step
    attr(cross$geno[[i]]$argmax,"off.end") <- off.end
  }

  cross
}

  

# end of argmax.geno.R
