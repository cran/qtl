#####################################################################
#
# util.R
#
# copyright (c) 2001-2, Karl W Broman, Johns Hopkins University
# last modified June, 2002
# first written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: pull.map, replace.map, create.map,
#           convert.cross, clean, drop.nullmarkers
#           drop.markers, geno.table, mf.k, mf.h, imf.k, imf.h
#           mf.cf, imf.cf, mf.m, imf.m, convert2ss, switch.order
#           subset.cross, fill.geno, check.covar
#           pull.chr (deprecated), find.marker, adjust.rf.ri,
#           getsex, getgenonames, fixXdata
#
######################################################################

######################################################################
#
# pull.map
#
# pull out the map portion of a cross object, as a list
#
######################################################################

pull.map <-
function(cross)
{
  a <- lapply(cross$geno,function(a) { b <- a$map; class(b) <- class(a); b })
  class(a) <- "map"
  a
  
}

######################################################################
#
# replace.map
#
# replace the map portion of a cross object with a list defining a map
#
######################################################################

replace.map <-
function(cross, map)
{
  n.chr <- nchr(cross) 
  n.mar <- nmar(cross)

  n.chr2 <- length(map)
  n.mar2 <- sapply(map,length)

  type <- class(cross)[1]
  if(type=="4way" || type=="f2ss") {
    mnames <- unlist(lapply(cross$geno, function(a) colnames(a$map)))
    mnames2 <- unlist(lapply(map, function(a) colnames(a)))
    n.mar2 <- n.mar2/2
  }
  else if(type == "bc" || type == "f2" || type == "riself" || type=="risib") {
    mnames <- unlist(lapply(cross$geno, function(a) names(a$map)))
    mnames2 <- unlist(lapply(map, function(a) names(a)))
  }
  else {
    stop(paste("Cross type", type, "not yet supported."))
  }

  # check that things line up properly
  if(n.chr != n.chr2)
    stop("Numbers of chromosomes don't match.")
  if(any(names(cross$geno) != names(map)))
    stop("Chromosome names don't match.")
  if(any(n.mar != n.mar2))
    stop("Number of markers don't match.")
  if(any(mnames != mnames2))
    stop("Marker names don't match.")

  # proceed if no errors
  for(i in 1:length(cross$geno))
    cross$geno[[i]]$map <- map[[i]]

  cross
}

######################################################################
#
# create.map
#
# create a new map with inserted inter-marker locations
#
# Note: map is a vector or a matrix with 2 rows
# 
######################################################################

create.map <-
function(map, step, off.end)
{
  if(step<0 || off.end<0) stop("step and off.end must be > 0.")

  if(!is.matrix(map)) { # sex-ave map
    if(length(map) == 1) { # just one marker!
      if(off.end==0) {
        if(step == 0) step <- 1
        nam <- names(map)
        map <- c(map,map+step)
        names(map) <- c(nam,paste("loc",step,sep=""))
      }
      else {
        if(step==0) m <- c(-off.end,off.end)
        else m <- seq(-off.end,off.end,by=step)
        m <- m[m!=0]
        names(m) <- paste("loc",m,sep="")
        map <- sort(c(m+map,map))
      }
      return(map)
    }

    minloc <- min(map)
    map <- map-minloc

    if(step==0 && off.end==0) return(map+minloc)
    else if(step==0 && off.end > 0) {
      a <- c(floor(min(map)-off.end),ceiling(max(map)+off.end))
      names(a) <- paste("loc", a, sep="")
      return(sort(c(a,map))+minloc)
    }
    else if(step>0 && off.end == 0) {
      a <- seq(floor(min(map)),max(map),
               by = step)
      if(any(is.na(match(a,map)))) {
        a <- a[is.na(match(a,map))]
        names(a) <- paste("loc",a,sep="")
        return(sort(c(a,map))+minloc)
      }
      else return(map+minloc)
    }
    else {
      a <- seq(floor(min(map)-off.end),ceiling(max(map)+off.end+step),
               by = step)
      a <- a[is.na(match(a,map))]
      
      # no more than one point above max(map)+off.end
      z <- (seq(along=a))[a >= max(map)+off.end]
      if(length(z) > 1) a <- a[-z[-1]]
      
      names(a) <- paste("loc",a,sep="")
      return(sort(c(a,map))+minloc)
    }
  } # end sex-ave map
  else { # sex-specific map
    minloc <- c(min(map[1,]),min(map[2,]))
    map <- map-minloc
    markernames <- colnames(map)

    if(step==0 && off.end==0) return(map+minloc)
    else if(step==0 && off.end > 0) {
      if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
        L1 <- L2 <- 1
      }
      else {
        L1 <- diff(range(map[1,]))
        L2 <- diff(range(map[2,]))
      }

      a <- c(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end))
      names(a) <- paste("loc", a, sep="")
      b <- c(floor(min(map[2,])-off.end)*L2/L1,
             ceiling(max(map[2,])+off.end)*L2/L1)
      n <- c(names(a)[1],markernames,names(a)[2])
      map <- cbind(c(a[1],b[1]),map,c(a[2],b[2]))
      dimnames(map) <- list(NULL,n)
      return(map+minloc)
    }
    else if(step>0 && off.end == 0) {
      if(ncol(map)==1) return(map+minloc)

      a <- seq(floor(min(map[1,])),max(map[1,]),
               by = step)
      a <- a[is.na(match(a,map[1,]))]
      
      b <- sapply(a,function(x,y,z) {
          ZZ <- min((seq(along=y))[y > x])
          (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1] }, map[1,],map[2,])
      m1 <- c(a,map[1,])
      m2 <- c(b,map[2,])
      names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
      return(rbind(sort(m1),sort(m2))+minloc)
    }
    else {
      a <- seq(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end+step),
               by = step)
      a <- a[is.na(match(a,map[1,]))]
      # no more than one point above max(map)+off.end
      z <- (seq(along=a))[a >= max(map[1,])+off.end]
      if(length(z) > 1) a <- a[-z[-1]]

      b <- sapply(a,function(x,y,z,ml) {
        if(x < min(y)) {
          return(min(z) - (min(y)-x)/diff(range(y))*diff(range(z)) - ml)
        }
        else if(x > max(y)) {
          return(max(z) + (x - max(y))/diff(range(y))*diff(range(z)) - ml)
        }
        else {
          ZZ <- min((seq(along=y))[y > x])
          (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1]
        }
        }, map[1,],map[2,], minloc[2])
      m1 <- c(a,map[1,])
      m2 <- c(b,map[2,])
      names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
      return(rbind(sort(m1),sort(m2))+minloc)
    }
  }
}

  
######################################################################
#
# convert.cross: convert a "qtl.cross" data set from the format
#                used in old versions (<= 0.65) of R/qtl to the
#                updated data structure (versions >= 0.70).
#
######################################################################

convert.cross <-
function(cross)   
{
  require(qtl)
  nchr <- length(cross$map)
  geno <- vector("list",nchr)
  nmar <- c(0,cumsum(sapply(cross$map,length)))
  for(i in 1:nchr) {
    whichpos <- (nmar[i]+1):nmar[i+1]
    geno[[i]] <- list(data=cross$geno[,whichpos],map=cross$map[[i]])
    dimnames(geno[[i]]$data) <- list(NULL, names(cross$map[[i]]))
    class(geno[[i]]) <- "A"
    chr.name <- names(cross$map)[i]
    if(chr.name == "X" || chr.name == "x")
      class(geno[[i]]) <- "X"
  }
  names(geno) <- names(cross$map)
  cross$geno <- geno
  type <- cross$type
  cross <- cross[1:2]
  class(cross) <- c(type,"cross")
  cross
}

    

######################################################################
#
# clean
# 
# remove all of the extraneous stuff from a cross object, to get back
# to just the data
#
######################################################################

clean <-
function(cross)
{
  cross2 <- list(geno=cross$geno,pheno=cross$pheno)

  for(i in 1:length(cross$geno)) {
    cross2$geno[[i]] <- list(data=cross$geno[[i]]$data,
                             map=cross$geno[[i]]$map)
    class(cross2$geno[[i]]) <- class(cross$geno[[i]])
  }
    
  class(cross2) <- class(cross)
  cross2
}


######################################################################
#
# drop.qtlgeno
#
# remove any QTLs from the genotype data and the genetic maps
# from data simulated via sim.cross. (They all have names "QTL*")
#
######################################################################

#drop.qtlgeno <-
#function(cross)  
#{
#  n.chr <- nchr(cross)
#  mar.names <- lapply(cross$geno, function(a) {
#    m <- a$map
#    if(is.matrix(m)) return(colnames(m))
#    else return(names(m)) } )
#    
#  for(i in 1:n.chr) {
#    o <- grep("^QTL[0-9]+",mar.names[[i]])
#    if(length(o) != 0) {
#      cross$geno[[i]]$data <- cross$geno[[i]]$data[,-o,drop=FALSE]
#      if(is.matrix(cross$geno[[i]]$map)) 
#        cross$geno[[i]]$map <- cross$geno[[i]]$map[,-o,drop=FALSE]
#      else
#        cross$geno[[i]]$map <- cross$geno[[i]]$map[-o]
#    }
#  }
#  cross
#}

######################################################################
#
# drop.nullmarkers
#
# remove markers that have no genotype data from the data matrix and
# genetic maps
#
######################################################################

drop.nullmarkers <-
function(cross)
{
  n.chr <- nchr(cross)

  keep.chr <- rep(TRUE,n.chr)
  for(i in 1:n.chr) {
    o <- !apply(cross$geno[[i]]$data,2,function(a) sum(!is.na(a)))
    if(any(o)) { # remove from genotype data and map
      mn.drop <- colnames(cross$geno[[i]]$data)[o]
      if(length(mn.drop) == ncol(cross$geno[[i]]$data)) 
        keep.chr[i] <- FALSE # removing all markers from this chromosome

      cross$geno[[i]]$data <- cross$geno[[i]]$data[,!o,drop=FALSE]

      if(is.matrix(cross$geno[[i]]$map)) 
        cross$geno[[i]]$map <- cross$geno[[i]]$map[,!o,drop=FALSE]
      else
        cross$geno[[i]]$map <- cross$geno[[i]]$map[!o]

      # results of calc.genoprob
      if(!is.na(match("prob",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$prob))
        cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,-o,,drop=FALSE]
      }

      # results of argmax.geno
      if(!is.na(match("argmax",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$argmax))
        cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-o,drop=FALSE]
      }

      # results of sim.geno
      if(!is.na(match("draws",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$draws))
        cross$geno[[i]]$draws <- cross$geno[[i]]$draws[,-o,,drop=FALSE]
      }

      # results of est.rf
      if(!is.na(match("rf",names(cross)))) {
        o <- match(mn.drop,colnames(cross$rf))
        cross$rf <- cross$rf[-o,-o]
      }
    }
  }
  cross$geno <- cross$geno[keep.chr]

  cross
}

    
######################################################################
#
# drop.markers
#
# remove a vector of markers from the data matrix and genetic maps
#
######################################################################

drop.markers <-
function(cross,markers)
{
  n.chr <- nchr(cross)

  keep.chr <- rep(TRUE,n.chr)
  for(i in 1:n.chr) {
    # find markers on this chromosome
    o <- match(markers,colnames(cross$geno[[i]]$data))
    o <- o[!is.na(o)]
    a <- rep(FALSE,ncol(cross$geno[[i]]$data))
    a[o] <- TRUE
    o <- a
    
    if(any(o)) { # remove from genotype data and map
      mn.drop <- colnames(cross$geno[[i]]$data)[o]
      if(length(mn.drop) == ncol(cross$geno[[i]]$data)) 
        keep.chr[i] <- FALSE # removing all markers from this chromosome

      cross$geno[[i]]$data <- cross$geno[[i]]$data[,!o,drop=FALSE]

      if(is.matrix(cross$geno[[i]]$map)) 
          x <- cross$geno[[i]]$map[,!o,drop=FALSE]
      else 
        cross$geno[[i]]$map <- cross$geno[[i]]$map[!o]

      # results of calc.genoprob
      if(!is.na(match("prob",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$prob))
        cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,-o,,drop=FALSE]
      }

      # results of argmax.geno
      if(!is.na(match("argmax",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$argmax))
        cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-o,drop=FALSE]
      }

      # results of sim.geno
      if(!is.na(match("draws",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$draws))
        cross$geno[[i]]$draws <- cross$geno[[i]]$draws[,-o,,drop=FALSE]
      }

      # results of est.rf
      if(!is.na(match("rf",names(cross)))) {
        o <- match(mn.drop,colnames(cross$rf))
        cross$rf <- cross$rf[-o,-o]
      }
    }
  }

  cross$geno <- cross$geno[keep.chr]

  cross
}

    
######################################################################
#
# geno.table
#
# create table showing observed numbers of individuals with each
# genotype at each marker
#
######################################################################

geno.table <- 
function(cross)
{
  n.chr <- nchr(cross)

  type <- class(cross)[1]
  if(type == "f2" || type=="f2ss") {
    n.gen <- 5
    gen.names <- c("AA","AB","BB","AA/AB","AB/BB")
  }
  else if(type == "bc") {
    n.gen <- 2
    gen.names <- c("AA","AB")
  }
  else if(type == "risib" || type=="riself") {
    n.gen <- 2
    gen.names <- c("AA","BB")
  }
  else if(type == "4way") {
    n.gen <- 10
    gen.names <- c("AC","BC","AD","BD","AC/AD","BC/BD",
                   "AC/BC","AD/BD","AC/BD","AD/BC")
  }
  else
    stop(paste("Unknown cross type:",type))
    
  res <- lapply(cross$geno, function(a,ngen) {
                a <- a$data; a[is.na(a)] <- 0
                apply(a,2,function(b,ngen) table(factor(b,levels=0:ngen)),ngen)
                },n.gen)
  results <- NULL
  for(i in 1:length(res)) 
    results <- rbind(results,t(res[[i]]))
  colnames(results) <- c("NA",gen.names)
  rownames(results) <- unlist(lapply(cross$geno,function(a) colnames(a$data)))
  results
}
    
# map functions
mf.k <- function(d) 0.5*tanh(d/50)
mf.h <- function(d) 0.5*(1-exp(-d/50))
imf.k <- function(r) 50*atanh(2*r)
imf.h <- function(r) -50*log(1-2*r)
mf.m <- function(d) sapply(d,function(a) min(a/100,0.5))
imf.m <- function(r) sapply(r,function(a) min(a*100,50))

# carter-falconer: mf.cf, imf.cf
imf.cf <- function(r) 12.5*(log(1+2*r)-log(1-2*r))+25*atan(2*r)

mf.cf <-
function(d)
{
  icf <- function(r,d)
    imf.cf(r)-d

  sapply(d,function(a) {
    if(a==0) return(0)
    uniroot(icf, c(0,0.5-1e-10),d=a)$root })
}


# convert F2 intercross to sex-specific
convert2ss <-
function(cross)  
{
  if(class(cross)[1] != "f2")
    stop("This function applies only to f2 crosses.")

  class(cross)[1] <- "f2ss"

  for(i in 1:nchr(cross)) 
    cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,
                                 cross$geno[[i]]$map)

  cross
}

######################################################################
#
# switch.order: change the marker order on a given chromosome to some
#               specified order
#
######################################################################

switch.order <-
function(cross, chr, order)
{
  # check chromosome argument
  if(is.character(chr)) {
    old.chr <- chr
    chr <- match(chr, names(cross$geno))
    if(length(chr) > 1) chr <- chr[1]
    if(is.na(chr))
      stop("There is no chromosome named", chr)
  }

  # check order argument
  n.mar <- nmar(cross)
  if(n.mar[chr] == length(order)-2 || n.mar[chr]==length(order)-1) 
    order <- order[1:n.mar[chr]]     # useful for output from ripple()
  if(n.mar[chr] != length(order))
    stop("Incorrect number of markers.")

  # remove any intermediate calculations, as they
  #   will no longer be meaningful
  cross <- clean(cross)

  # re-order markers
  cross$geno[[chr]]$data <- cross$geno[[chr]]$data[,order,drop=FALSE]
  m <- seq(0,by=5,length=ncol(cross$geno[[chr]]$data))
  names(m) <- colnames(cross$geno[[chr]]$data)
  if(is.matrix(cross$geno[[chr]]$map)) 
    cross$geno[[chr]]$map <- rbind(m,m)
  else
    cross$geno[[chr]]$map <- m

  # re-estimate map
  newmap <- est.map(subset(cross,chr=chr))
  cross$geno[[chr]]$map <- newmap[[1]]

  cross
}

######################################################################
#
# subset.cross: General subsetting function for a cross object
#
######################################################################

subset.cross <-
function(x, chr, ind, ...)  
{
  n.chr <- nchr(x)
  n.ind <- nind(x)

  # pull out relevant chromosomes
  if(!missing(chr)) {
    if(is.logical(chr)) {
      if(length(chr) != n.chr)
        stop(paste("If logical, chr argument must have length", n.chr))
      chr <- (1:n.chr)[chr]
    }
        
    if(is.numeric(chr)) {
      # if all negative numbers, convert to positive
      if(all(chr < 1)) chr <- (1:n.chr)[chr]
        
      if(any(chr < 1 | chr > n.chr))
        stop("Chromosome numbers out of range.")
    }
    else {
      if(any(is.na(match(chr,names(x$geno)))))
        stop("Not all chromosome names found.")
      # convert to numeric
      chr <- match(chr,names(x$geno))
    }

    if(!is.na(match("rf",names(x)))) { # pull out part of rec fracs
      n.mar <- nmar(x)
      n.chr <- n.chr
      wh <- rbind(c(0,cumsum(n.mar)[-n.chr])+1,cumsum(n.mar))
      dimnames(wh) <- list(NULL, names(n.mar))
      wh <- wh[,chr,drop=FALSE]
      wh <- unlist(apply(wh,2,function(a) a[1]:a[2]))
      x$rf <- x$rf[wh,wh]
    }

    x$geno <- x$geno[chr]
  }

  if(!missing(ind)) {
    if(is.logical(ind)) {
      if(length(ind) != n.ind)
        stop(paste("If logical, ind argument must have length", n.ind))
      ind <- (1:n.ind)[ind]
    }
        
    if(is.numeric(ind)) {
      # if all negative numbers, convert to positive
      if(all(ind < 1)) ind <- (1:n.ind)[ind]
        
      if(any(ind < 1 | ind > n.ind))
        stop("Individual numbers out of range.")
    }
    else stop("ind argument must be either logical or numeric.")
    # Note: ind should now be a numeric vector

    if(length(ind) == 0)
      stop("Must retain at least one individual.")
    if(length(ind) == 1)
      warning("Retained only one individual!")

    x$pheno <- x$pheno[ind,,drop=FALSE]

    for(i in 1:nchr(x)) {
      x$geno[[i]]$data <- x$geno[[i]]$data[ind,,drop=FALSE]

      if(!is.na(match("prob",names(x$geno[[i]])))) 
        x$geno[[i]]$prob <- x$geno[[i]]$prob[ind,,,drop=FALSE]
      if(!is.na(match("errorlod",names(x$geno[[i]])))) 
        x$geno[[i]]$errorlod <- x$geno[[i]]$errorlod[ind,,drop=FALSE]
      if(!is.na(match("argmax",names(x$geno[[i]])))) 
        x$geno[[i]]$argmax <- x$geno[[i]]$argmax[ind,,drop=FALSE]
      if(!is.na(match("draws",names(x$geno[[i]])))) 
        x$geno[[i]]$draws <- x$geno[[i]]$draws[ind,,,drop=FALSE]
    }
  }
  x
}

pull.chr <-
function(cross, chr) {
  warning("pull.chr is deprecated; use subset.cross.")
  subset.cross(cross, chr)
}


######################################################################
#
# c.cross: Combine crosses
#
######################################################################

c.cross <-
function(...)
{
  args <- list(...)

  x <- args[[1]]
  if(class(x)[2] != "cross")
    stop("All arguments must be cross arguments.")
  n.phe <- nphe(x)
  phenam <- colnames(x$pheno)
  n.chr <- nchr(x)
  chrnam <- names(x$geno)
  n.mar <- nmar(x)
  marnam <- unlist(lapply(x$geno,function(b) names(b$map)))
  type <- class(x)[1]
  x <- clean(x)

  # if only one cross, just return it
  if(length(args)==1) return(args[[1]])

  for(i in 2:length(args)) {
    y <- args[[i]]
    if(class(y)[2] != "cross")
      stop("All arguments must be cross arguments.")
    if(type != class(y)[1])
      stop("All arguments must be the same cross type.")
    if(n.phe != nphe(y) || any(phenam != colnames(y$pheno))) 
      stop("All arguments must have the same phenotypes.")
    if(n.chr != nchr(y) || any(chrnam != names(y$geno)) ||
       any(n.mar != nmar(y)) ||
       any(marnam != unlist(lapply(y$geno,function(b) names(b$map)))))
      stop("All arguments have have the same chromosomes and markers.")

    x$pheno <- rbind(x$pheno,y$pheno)
    for(j in 1:n.chr) 
      x$geno[[j]]$data <- rbind(x$geno[[j]]$data,y$geno[[j]]$data)
  }

  x
}

######################################################################
#
# fill.geno: Run argmax.geno or sim.geno and then fill in the
#            genotype data with the results.  This will allow
#            rough genome scans by marker regression without
#            holes.  WE WOULD NOT PLACE ANY TRUST IN THE RESULTS!
#
######################################################################

fill.geno <-
function(cross, method=c("imp","argmax"), error.prob=0,
         map.function=c("haldane","kosambi","c-f","morgan"))
{
  method <- match.arg(method)
  
  # don't let error.prob be exactly zero (or >1)
  if(error.prob < 1e-50) error.prob <- 1e-50
  if(error.prob > 1) {
    error.prob <- 1-1e-50
    warning("error.prob shouldn't be > 1!")
  }

  # remove any extraneous material
  cross <- clean(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)

  if(method=="imp") {
    # do one imputation
    temp <- sim.geno(cross,n.draws=1,step=0,off.end=0,
                     error.prob=error.prob,map.function=map.function)
    # replace the genotype data with the results,
    #     stripping off any attributes
    for(i in 1:n.chr) {
      nam <- colnames(cross$geno[[i]]$data)
      cross$geno[[i]]$data <-
        matrix(as.numeric(temp$geno[[i]]$draws[,,1]),ncol=n.mar[i])
      colnames(cross$geno[[i]]$data) <- nam
    }
  }
  else {
    # run the Viterbi algorithm
    temp <- argmax.geno(cross,step=0,off.end=0,error.prob=error.prob,
                        map.function=map.function)
    # replace the genotype data with the results,
    #     stripping off any attributes
    for(i in 1:n.chr) {
      nam <- colnames(cross$geno[[i]]$data)
      cross$geno[[i]]$data <-
        matrix(as.numeric(temp$geno[[i]]$argmax),ncol=n.mar[i])
      colnames(cross$geno[[i]]$data) <- nam
    }
  }
  cross
}

######################################################################
#
# checkcovar
#
# This is a utility function for scanone and scantwo.  We remove  
# individuals with missing phenotypes or covariates and check
# that the covariates are of the appropriate form.
#
######################################################################

checkcovar <-
function(cross, pheno.col, addcovar, intcovar)
{
  # check phenotypes
  if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
  if(pheno.col < 1 || pheno.col > nphe(cross))
    stop("Specified phenotype column is invalid.")

  orig.n.ind <- nind(cross)

  # drop individuals with missing phenotypes
  pheno <- cross$pheno[,pheno.col]
  if(any(is.na(pheno))) {
    keep.ind <- (1:length(pheno))[!is.na(pheno)]
    cross <- subset.cross(cross,ind=keep.ind)
    pheno <- pheno[keep.ind]
  }
  else keep.ind <- 1:nind(cross)
  n.ind <- nind(cross)
  n.chr <- nchr(cross)      # number of chromosomes
  type <- class(cross)[1]   # type of cross

  n.addcovar <- n.intcovar <- 0
  if(!is.null(addcovar)) { # for additive covariates
    if(!is.matrix(addcovar)) {
      if(is.vector(addcovar) || is.data.frame(addcovar))
        addcovar <- as.matrix(addcovar)
      else stop("addcovar should be a matrix")
    }
    if(!all(apply(addcovar,2,is.numeric)))
      stop("All columns of addcovar must be numeric")
    if( nrow(addcovar) != orig.n.ind ) {
      # the length of additive covariates is incorrect
      stop("Number of rows in additive covariates is incorrect")
    }
    addcovar <- addcovar[keep.ind,,drop=FALSE]
    n.addcovar <- ncol(addcovar)
  }
  if(!is.null(intcovar)) { # interacting covariates
    if(!is.matrix(intcovar)) {
      if(is.vector(intcovar) || is.data.frame(intcovar))
        intcovar <- as.matrix(intcovar)
      else stop("intcovar should be a matrix")
    }
    if(!all(apply(intcovar,2,is.numeric)))
      stop("All columns of intcovar must be numeric")
    if(nrow(intcovar)[1] != orig.n.ind) {
      # the length of interacting covariates is incorrect
      stop("The length of interacting covariates is incorrect!")
    }
    intcovar <- intcovar[keep.ind,,drop=FALSE]
    n.intcovar <- ncol(intcovar)
  }

  # drop individuals missing any covariates
  if(!is.null(addcovar)) { # note that intcovar is contained in addcovar
    wh <- apply(cbind(addcovar,intcovar),1,function(a) any(is.na(a)))
    if(any(wh)) {
      cross <- subset.cross(cross,ind=(!wh))
      pheno <- pheno[!wh]
      addcovar <- addcovar[!wh,,drop=FALSE]
      if(!is.null(intcovar)) intcovar <- intcovar[!wh,,drop=FALSE]
      n.ind <- nind(cross)
      warning("Dropping individuals with missing covariates")
    }
  }

  # make sure columns of intcovar are contained in addcovar
  if(!is.null(intcovar)) {
    if(is.null(addcovar)) {
      addcovar <- intcovar
      n.addcovar <- n.intcovar
      warning("addcovar forced to contain all columns of intcovar")
    }
    else {
      wh <- 1:n.intcovar
      for(i in 1:n.intcovar) {
        o <- (apply(addcovar,2,function(a,b) max(abs(a-b)),intcovar[,i])<1e-14)
        if(any(o)) wh[i] <- (1:n.addcovar)[o]
        else wh[i] <- NA
      }
      if(any(is.na(wh))) {
        addcovar <- cbind(addcovar,intcovar[,is.na(wh)])
        n.addcovar <- ncol(addcovar)
        warning("addcovar forced to contain all columns of intcovar")
      }
    }
  }

  list(cross=cross, pheno=pheno, addcovar=addcovar, intcovar=intcovar,
       n.addcovar=n.addcovar, n.intcovar=n.intcovar)
}

# Find the nearest marker to a particular position
find.marker <-
function(cross, chr, pos)  
{
  # if chr has length 1, expand if necessary
  if(length(chr) == 1) 
    chr <- rep(chr,length(pos))
  # otherwise, chr and pos should have same length
  else if(length(chr) != length(pos)) 
    stop("chr and pos must be the same length.")

  markers <- rep("",length(chr))
  for(i in 1:length(chr)) {
    # find chromosome
    o <- match(chr[i], names(cross$geno))
    if(is.na(o)) markers[i] <- NA  # chr not matched
    else {
      thismap <- cross$geno[[o]]$map # genetic map

      # sex-specific map; look at female positions
      if(is.matrix(thismap)) thismap <- thismap[1,]
      
      # find closest marker
      d <- abs(thismap-pos[i])
      o2 <- (1:length(d))[d==min(d)]
      if(length(o2)==1) markers[i] <- names(thismap)[o2]
      # if multiple markers are equidistant,
      #     choose the one with the most data
      #     or choose among them at random
      else {
        x <- names(thismap)[o2]
        n.geno <- apply(cross$geno[[o]]$data[,o2],2,function(a) sum(!is.na(a)))
        o2 <- o2[n.geno==max(n.geno)]
        markers[i] <- names(thismap)[sample(o2,1)]
      }
    }
  }

  markers
}


# expand recombination fractions for RI lines
adjust.rf.ri <-
function(r, type=c("self","sib"), expand=TRUE)
{
  # type of RI lines
  type <- match.arg(type)

  if(type=="self") {
    if(expand) return(r*2/(1+2*r))
    else return(r/2/(1-r))
  }
  else {
    if(expand) return(r*4/(1+6*r))
    else return(r/(4-6*r))
  }
}

# get sex and pgm columns from phenotype data
getsex <-
function(cross)
{
  phe.names <- names(cross$pheno)

  sex.column <- grep("^[Ss][Ee][Xx]$", phe.names)
  pgm.column <- grep("^[Pp][Gg][Mm]$", phe.names)

  if(length(sex.column)==0) { # no sex included
    sex <- NULL
  }
  else {
    temp <- cross$pheno[,sex.column[1]]
    if(is.numeric(temp)) {
      if(any(!is.na(temp) & temp != 0 & temp != 1)) {
        warning("Sex column should be coded as 0=female 1=male; sex ignored.")
        sex <= NULL
      }
      else sex <- temp
    }
    else {
      if(!is.factor(temp)) temp <- as.factor(temp)

      if(length(levels(temp)) != 2) {
        warning("Sex column should be coded as a two-level factor; sex ignored.")
        sex <- NULL
      }
      else { # is a factor with two levels
        lev <- levels(temp)
        if(length(grep("^[Ff]",lev))>0 &&
           length(males <- grep("^[Mm]",lev))>0) {
          temp <- as.character(temp)
          sex <- rep(0,length(temp))
          sex[is.na(temp)] <- NA
          sex[!is.na(temp) & temp==lev[males]] <- 1
        }
        else 
          warning("Don't understand levels in sex column; sex ignored.")
      }
    }
  }

  if(length(pgm.column)==0) { # no sex included
    pgm <- NULL
  }
  else {
    temp <- cross$pheno[,pgm.column[1]]
    if(!is.numeric(temp))
      temp <- as.numeric(as.character(pgm))
    if(any(!is.na(temp) & temp != 0 & temp != 1)) {
      warning("pgm column should be coded as 0/1; pgm ignored.")
      pgm <- NULL
    }
    else pgm <- temp
  }

  list(sex=sex,pgm=pgm)
}
          


# get names of genotypes
getgenonames <-
function(type=c("f2","bc","f2ss","riself","risib","4way"),
         chrtype=c("A","X"), x.treatment=c("simple","full"),
         sexpgm)
{  
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  x.treatment <- match.arg(x.treatment)

  if(type == "f2" || type=="f2ss") {
    if(chrtype == "A")  # autosomal
      gen.names <- c("AA","AB","BB")
    else { # X chromsome
      if(length(sexpgm$sex)>0) {
        if(length(sexpgm$pgm)>0) {
          if(x.treatment=="simple") 
            gen.names <- c("AA","AB","BB")
          else {
            gen.names <- c("AA","AB","BB","AY","BY")
          }
        }
        else { # no pgm
          if(x.treatment=="simple") 
            gen.names <- c("AA","AB","BB")
          else 
            gen.names <- c("AA","AB","AY","BY")
        }
      }
      else { # no sex
        if(length(sexpgm$pgm)>0) 
          gen.names <- c("AA","AB","BB")
        else
          gen.names <- c("AA","AB")
      }
    }
  }
  else if(type == "bc") {
    if(chrtype=="X" && length(sexpgm$sex)>0) {
      if(x.treatment=="simple") 
        gen.names <- c("AA","AB","BB")
      else
        gen.names <- c("AA","AB","AY","BY")
    }
    else gen.names <- c("AA","AB")
  }
  else if(type == "riself" || type=="risib") 
    gen.names <- c("AA","BB")
  else if(type == "4way") 
    gen.names <- c("AC","BC","AD","BD")
  else stop(paste("Cross type", type, "not supported."))

  gen.names
}

# revise genotype data, probabilities or imputations for the X chromosome
fixXdata <-
function(type=c("f2ss","f2","bc"), x.treatment=c("simple","full"),
         sexpgm, geno, prob, draws, pairprob)
{
  type <- match.arg(type)
  x.treatment <- match.arg(x.treatment)

  sex <- sexpgm$sex
  pgm <- sexpgm$pgm

  if(!missing(geno)) {
    if(type == "f2" || type=="f2ss") {
      if(length(sex)>0) {
        if(length(pgm)>0) {
          if(x.treatment=="simple") {
            w <- geno[sex==1,,drop=FALSE]
            w[!is.na(w) & w==2] <- 3
            geno[sex==1,] <- w
            
            w <- geno[sex==0 & pgm==1,,drop=FALSE]
            w[!is.na(w) & w==1] <- 3
            geno[sex==0 & pgm==1,] <- w
          }
          else  {
            w <- geno[sex==1,,drop=FALSE]
            w[!is.na(w) & w==1] <- 4
            w[!is.na(w) & w==2] <- 5
            geno[sex==1,] <- w
            
            w <- geno[sex==0 & pgm==1,,drop=FALSE]
            w[!is.na(w) & w==1] <- 3
            geno[sex==0 & pgm==1,] <- w
          }
        }
        else { # no pgm
          if(x.treatment=="simple") {
            w <- geno[sex==1,,drop=FALSE]
            w[!is.na(w) & w==2] <- 3
            geno[sex==1,] <- w
          }
          else {
            w <- geno[sex==1,,drop=FALSE]
            w[!is.na(w) & w==1] <- 3
            w[!is.na(w) & w==2] <- 4
            geno[sex==1,] <- w
          }
        }
      }
      else { # no sex; assume all females
        if(length(pgm)>0) { 
          w <- geno[pgm==1,,drop=FALSE]
          w[!is.na(w) & w==1] <- 3
          geno[pgm==1,] <- w
        }
        # else leave unchanged
      }
    } # end f2
    else { # backcross
      if(length(sex)>0) {
        if(x.treatment=="simple") {
          w <- geno[sex==1,,drop=FALSE]
          w[!is.na(w) & w==2] <- 3
          geno[sex==1,] <- w
        }
        else {
          w <- geno[sex==1,,drop=FALSE]
          w[!is.na(w) & w==1] <- 3
          w[!is.na(w) & w==2] <- 4
          geno[sex==1,] <- w
        }
      }
      # else leave unchanged
    } # end backcross
    return(geno)
  } # end if(!missing(geno))

  else if(!missing(prob)) {
    n.gen <- length(getgenonames(type,"X",x.treatment,sexpgm))
    newdim <- dim(prob)
    newdim[3] <- n.gen
    newprob <- array(0,dim=newdim)
    
    if(type == "f2" || type=="f2ss") {
      if(length(sex)>0) {
        if(length(pgm)>0) {
          if(x.treatment=="simple") {
            newprob[sex==1,,c(1,3)] <- prob[sex==1,,]
            newprob[sex==0 & pgm==0,,1:2] <- prob[sex==0 & pgm==0,,]
            newprob[sex==0 & pgm==1,,c(3,2)] <- prob[sex==0 & pgm==1,,]
          }
          else  {
            newprob[sex==1,,4:5] <- prob[sex==1,,]
            newprob[sex==0 & pgm==0,,1:2] <- prob[sex==0 & pgm==0,,]
            newprob[sex==0 & pgm==1,,c(3,2)] <- prob[sex==0 & pgm==1,,]
          }
        }
        else { # no pgm
          if(x.treatment=="simple") {
            newprob[sex==1,,c(1,3)] <- prob[sex==1,,]
            newprob[sex==0,,1:2] <- prob[sex==0,,]
          }
          else {
            newprob[sex==1,,3:4] <- prob[sex==1,,]
            newprob[sex==0,,1:2] <- prob[sex==0,,]
          }
        }
      }
      else { # no sex; assume all females
        if(length(pgm)>0) { 
            newprob[pgm==0,,1:2] <- prob[pgm==0,,]
            newprob[pgm==1,,c(3,2)] <- prob[pgm==1,,]
        }
        else newprob <- prob
      }
    } # end f2
    else { # backcross
      if(length(sex)>0) {
        if(x.treatment=="simple") {
          newprob[sex==1,,c(1,3)] <- prob[sex==1,,]
          newprob[sex==0,,1:2] <- prob[sex==0,,]
        }
        else {
          newprob[sex==1,,3:4] <- prob[sex==1,,]
          newprob[sex==0,,1:2] <- prob[sex==0,,]
        }
      }
      else newprob <- prob
    } # end backcross
    
    attr(newprob,"error.prob") <- attr(prob,"error.prob")
    attr(newprob,"step") <- attr(prob,"step")
    attr(newprob,"off.end") <- attr(prob,"off.end")
    attr(newprob,"map.function") <- attr(prob,"map.function")
    return(newprob)
  }

  else if(!missing(draws)) {
    if(type == "f2" || type=="f2ss") {
      if(length(sex)>0) {
        if(length(pgm)>0) {
          if(x.treatment=="simple") {
            w <- draws[sex==1,,,drop=FALSE]
            w[!is.na(w) & w==2] <- 3
            draws[sex==1,,] <- w
            
            w <- draws[sex==0 & pgm==1,,,drop=FALSE]
            w[!is.na(w) & w==1] <- 3
            draws[sex==0 & pgm==1,,] <- w
          }
          else  {
            w <- draws[sex==1,,,drop=FALSE]
            w[!is.na(w) & w==1] <- 4
            w[!is.na(w) & w==2] <- 5
            draws[sex==1,,] <- w
            
            w <- draws[sex==0 & pgm==1,,,drop=FALSE]
            w[!is.na(w) & w==1] <- 3
            draws[sex==0 & pgm==1,,] <- w
          }
        }
        else { # no pgm
          if(x.treatment=="simple") {
            w <- draws[sex==1,,,drop=FALSE]
            w[!is.na(w) & w==2] <- 3
            draws[sex==1,,] <- w
          }
          else {
            w <- draws[sex==1,,,drop=FALSE]
            w[!is.na(w) & w==1] <- 3
            w[!is.na(w) & w==2] <- 4
            draws[sex==1,,] <- w
          }
        }
      }
      else { # no sex; assume all females
        if(length(pgm)>0) { 
          w <- draws[pgm==1,,,drop=FALSE]
          w[!is.na(w) & w==1] <- 3
          draws[pgm==1,,] <- w
        }
        # else leave unchanged
      }
    } # end f2
    else { # backcross
      if(length(sex)>0) {
        if(x.treatment=="simple") {
          w <- draws[sex==1,,,drop=FALSE]
          w[!is.na(w) & w==2] <- 3
          draws[sex==1,,] <- w
        }
        else {
          w <- draws[sex==1,,,drop=FALSE]
          w[!is.na(w) & w==1] <- 3
          w[!is.na(w) & w==2] <- 4
          draws[sex==1,,] <- w
        }
      }
      # else leave unchanged
    } # end backcross
    return(draws)
  }

  else if(!missing(pairprob)) {
    n.gen <- length(getgenonames(type,"X",x.treatment,sexpgm))
    newdim <- dim(pairprob)
    newdim[3:4] <- n.gen
    newpairprob <- array(0,dim=newdim)
    
    if(type == "f2" || type=="f2ss") {
      if(length(sex)>0) {
        if(length(pgm)>0) {
          if(x.treatment=="simple") {
            newpairprob[sex==1,,c(1,3),c(1,3)] <- pairprob[sex==1,,,]
            newpairprob[sex==0 & pgm==0,,1:2,1:2] <- pairprob[sex==0 & pgm==0,,,]
            newpairprob[sex==0 & pgm==1,,c(3,2),c(3,2)] <- pairprob[sex==0 & pgm==1,,,]
          }
          else  {
            newpairprob[sex==1,,4:5,4:5] <- pairprob[sex==1,,,]
            newpairprob[sex==0 & pgm==0,,1:2,1:2] <- pairprob[sex==0 & pgm==0,,,]
            newpairprob[sex==0 & pgm==1,,c(3,2),c(3,2)] <- pairprob[sex==0 & pgm==1,,,]
          }
        }
        else { # no pgm
          if(x.treatment=="simple") {
            newpairprob[sex==1,,c(1,3),c(1,3)] <- pairprob[sex==1,,,]
            newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]
          }
          else {
            newpairprob[sex==1,,3:4,3:4] <- pairprob[sex==1,,,]
            newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]
          }
        }
      }
      else { # no sex; assume all females
        if(length(pgm)>0) { 
            newpairprob[pgm==0,,1:2,1:2] <- pairprob[pgm==0,,,]
            newpairprob[pgm==1,,c(3,2),c(3,2)] <- pairprob[pgm==1,,,]
        }
        else newpairprob <- pairprob
      }
    } # end f2
    else { # backcross
      if(length(sex)>0) {
        if(x.treatment=="simple") {
          newpairprob[sex==1,,c(1,3),c(1,3)] <- pairprob[sex==1,,,]
          newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]
        }
        else {
          newpairprob[sex==1,,3:4,3:4] <- pairprob[sex==1,,,]
          newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]
        }
      }
      else newpairprob <- pairprob
    } # end backcross
    
    return(newpairprob)
  }

  else stop("Need to specify geno, prob, draws, or pairprob")
}

# end of util.R

