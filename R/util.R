######################################################################
#
# util.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# Oct, 2001; Sept, 2001; July, 2001; Apr, 2001; Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: pull.map, replace.map, pull.chr, create.map,
#           convert.cross, clean, drop.qtlgeno, drop.nullmarkers
#           drop.markers, geno.table, mf.k, mf.h, imf.k, imf.h
#           mf.cf, imf.cf, convert2ss, switch.order
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
  else if(type == "bc" || type == "f2") {
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
# pull.chr
#
# Pull out a portion of the chromosomes from a cross object
#
######################################################################

pull.chr <-
function(cross, chr)
{
  if(is.numeric(chr)) {
    if(any(chr < 1 | chr > nchr(cross)))
      stop("Chromosome numbers out of range.")
  }
  else {
    if(any(is.na(match(chr,names(cross$geno)))))
      stop("Not all chromosome names found.")
  }

  if(!is.na(match("rf",names(cross)))) { # pull out part of rec fracs
    n.mar <- nmar(cross)
    n.chr <- nchr(cross)
    wh <- rbind(c(0,cumsum(n.mar)[-n.chr])+1,cumsum(n.mar))
    dimnames(wh) <- list(NULL, names(n.mar))
    wh <- as.matrix(wh[,chr])
    wh <- unlist(apply(wh,2,function(a) a[1]:a[2]))
    cross$rf <- cross$rf[wh,wh]
  }

  cross$geno <- cross$geno[chr]
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
  if(!is.matrix(map)) { # sex-ave map
    if(step==0 && off.end==0) return(map)
    else if(step==0 && off.end > 0) {
      a <- c(floor(min(map)-off.end),ceiling(max(map)+off.end))
      names(a) <- paste("loc", a, sep="")
      return(sort(c(a,map)))
    }
    else if(step>0 && off.end == 0) {
      a <- seq(floor(min(map)),max(map),
               by = step)
      if(any(is.na(match(a,map)))) {
        a <- a[is.na(match(a,map))]
        names(a) <- paste("loc",a,sep="")
        return(sort(c(a,map)))
      }
      else return(map)
    }
    else {
      a <- seq(floor(min(map)-off.end),ceiling(max(map)+off.end+step),
               by = step)
      a <- a[is.na(match(a,map))]
      
      # no more than one point above max(map)+off.end
      z <- (1:length(a))[a >= max(map)+off.end]
      if(length(z) > 1) a <- a[-z[-1]]
      
      names(a) <- paste("loc",a,sep="")
      return(sort(c(a,map)))
    }
  } # end sex-ave map
  else { # sex-specific map
    if(step==0 && off.end==0) return(map)
    else if(step==0 && off.end > 0) {
      L1 <- diff(range(map[1,]))
      L2 <- diff(range(map[2,]))
      a <- c(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end))
      names(a) <- paste("loc", a, sep="")
      b <- c(floor(min(map[2,])-off.end)*L2/L1,
             ceiling(max(map[1,])+off.end)*L2/L1)
      n <- c(names(a)[1],colnames(map),names(a)[2])
      map <- cbind(c(a[1],b[1]),map,c(a[2],b[2]))
      dimnames(map) <- list(NULL,n)
      return(map)
    }
    else if(step>0 && off.end == 0) {
      a <- seq(floor(min(map[1,])),max(map[1,]),
               by = step)
      a <- a[is.na(match(a,map[1,]))]
      names(a) <- paste("loc",a,sep="")
      b <- sapply(a,function(x,y,z) {
          I <- min((1:length(y))[y > x])
          (x-y[I-1])/(y[I]-y[I-1])*(z[I]-z[I-1])+z[I-1] }, map[1,],map[2,])
      names(b) <- names(a)
      return(rbind(sort(c(a,map[1,])),sort(c(b,map[2,]))))
    }
    else {
      a <- seq(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end+step),
               by = step)
      a <- a[is.na(match(a,map[1,]))]
      # no more than one point above max(map)+off.end
      z <- (1:length(a))[a >= max(map[1,])+off.end]
      if(length(z) > 1) a <- a[-z[-1]]
      names(a) <- paste("loc",a,sep="")

      b <- sapply(a,function(x,y,z) {
        if(x < min(y))
          return( min(z) - (min(y)-x)/diff(range(y))*diff(range(z)) )
        else if(x > max(y))
          return( max(z) + (x - max(y))/diff(range(y))*diff(range(z)) )
        else {
          I <- min((1:length(y))[y > x])
          (x-y[I-1])/(y[I]-y[I-1])*(z[I]-z[I-1])+z[I-1]
        }
        }, map[1,],map[2,])
      names(b) <- names(a)

      return(rbind(sort(c(a,map[1,])), sort(c(b,map[2,]))))
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

drop.qtlgeno <-
function(cross)  
{
  n.chr <- nchr(cross)
  mar.names <- lapply(cross$geno, function(a) {
    m <- a$map
    if(is.matrix(m)) return(colnames(m))
    else return(names(m)) } )
    
  for(i in 1:n.chr) {
    o <- grep("^QTL[0-9]+",mar.names[[i]])
    if(length(o) != 0) {
      cross$geno[[i]]$data <- cross$geno[[i]]$data[,-o]
      if(is.matrix(cross$geno[[i]]$map)) 
        cross$geno[[i]]$map <- cross$geno[[i]]$map[,-o]
      else
        cross$geno[[i]]$map <- cross$geno[[i]]$map[-o]
    }
  }
  cross
}

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

      if(sum(!o) == 1) mn <- colnames(cross$geno[[i]]$data)[!o]

      cross$geno[[i]]$data <- cross$geno[[i]]$data[,!o]

      if(is.matrix(cross$geno[[i]]$map)) {
        if(sum(!o) == 1) {
          x <- as.matrix(cross$geno[[i]]$map[,!o])
          colnames(x) <- mn
        }
        else 
          cross$geno[[i]]$map <- cross$geno[[i]]$map[,!o]
      }
      else 
        cross$geno[[i]]$map <- cross$geno[[i]]$map[!o]

      if(sum(!o) == 1) {
        cross$geno[[i]]$data <- as.matrix(cross$geno[[i]]$data)
        colnames(cross$geno[[i]]$data) <- mn
      }

      # results of calc.genoprob
      if(!is.na(match("prob",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$prob))
        cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,-o,]
      }

      # results of argmax.geno
      if(!is.na(match("argmax",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$argmax))
        cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-o]
      }

      # results of sim.geno
      if(!is.na(match("draws",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$draws))
        cross$geno[[i]]$draws <- cross$geno[[i]]$draws[,-o,]
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

      if(sum(!o) == 1) mn <- colnames(cross$geno[[i]]$data)[!o]

      cross$geno[[i]]$data <- cross$geno[[i]]$data[,!o]

      if(is.matrix(cross$geno[[i]]$map)) {
        if(sum(!o) == 1) {
          x <- as.matrix(cross$geno[[i]]$map[,!o])
          colnames(x) <- mn
        }
        else 
          cross$geno[[i]]$map <- cross$geno[[i]]$map[,!o]
      }
      else 
        cross$geno[[i]]$map <- cross$geno[[i]]$map[!o]

      if(sum(!o) == 1) {
        cross$geno[[i]]$data <- as.matrix(cross$geno[[i]]$data)
        colnames(cross$geno[[i]]$data) <- mn
      }

      # results of calc.genoprob
      if(!is.na(match("prob",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$prob))
        cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,-o,]
      }

      # results of argmax.geno
      if(!is.na(match("argmax",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$argmax))
        cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-o]
      }

      # results of sim.geno
      if(!is.na(match("draws",names(cross$geno[[i]])))) {
        o <- match(mn.drop,colnames(cross$geno[[i]]$draws))
        cross$geno[[i]]$draws <- cross$geno[[i]]$draws[,-o,]
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
  if(type == "f2") {
    n.gen <- 5
    gen.names <- c("AA","AB","BB","AA/AB","AB/BB")
  }
  else if(type == "bc") {
    n.gen <- 2
    gen.names <- c("AA","AB")
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
  if(n.mar[chr] == length(order)-2) # useful for output from ripple()
    order <- order[1:n.mar[chr]]
  if(n.mar[chr] != length(order))
    stop("Incorrect number of markers.")

  # remove any intermediate calculations, as they
  #   will no longer be meaningful
  cross <- clean(cross)

  # re-order markers
  cross$geno[[chr]]$data <- cross$geno[[chr]]$data[,order]
  m <- seq(0,by=5,length=ncol(cross$geno[[chr]]$data))
  names(m) <- colnames(cross$geno[[chr]]$data)
  if(is.matrix(cross$geno[[chr]]$map)) 
    cross$geno[[chr]]$map <- rbind(m,m)
  else
    cross$geno[[chr]]$map <- m

  # re-estimate map
  newmap <- est.map(pull.chr(cross,chr))
  cross$geno[[chr]]$map <- newmap[[1]]

  cross
}

# end of util.R
