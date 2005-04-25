#####################################################################
#
# util.R
#
# copyright (c) 2001-5, Karl W Broman, Johns Hopkins University
#     [find.pheno and find.flanking from Brian Yandell]
# last modified Apr, 2005
# first written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: pull.map, replace.map, create.map,
#           convert.cross, clean, drop.nullmarkers
#           drop.markers, geno.table, mf.k, mf.h, imf.k, imf.h
#           mf.cf, imf.cf, mf.m, imf.m, convert2ss, switch.order
#           subset.cross, fill.geno, check.covar, find.marker,
#           adjust.rf.ri, pull.geno, lodint, bayesint, makeSSmap,
#           comparecrosses, movemarker, summary.map,
#           print.summary.map, convert.scanone, find.pheno,
#           find.flanking, strip.partials, comparegeno
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
  else 
    stop("Cross type ", type, " not yet supported.")

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
  found <- rep(FALSE, length(markers))
  for(i in 1:n.chr) {
    # find markers on this chromosome
    o <- match(markers,colnames(cross$geno[[i]]$data))
    found[!is.na(o)] <- TRUE
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

  if(any(!found)) 
    warning("Markers not found: ", paste(markers[!found],collapse=" "))

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
  else stop("Unknown cross type: ",type)
    
  res <- lapply(cross$geno, function(a,ngen) {
                a <- a$data; a[is.na(a)] <- 0
                apply(a,2,function(b,ngen) table(factor(b,levels=0:ngen)),ngen)
                },n.gen)
  results <- NULL
  for(i in 1:length(res)) 
    results <- rbind(results,t(res[[i]]))
  colnames(results) <- c("missing",gen.names)
  rownames(results) <- unlist(lapply(cross$geno,function(a) colnames(a$data)))

  pval <- rep(NA,nrow(results))
  if(type=="bc" || type=="risib" || type=="riself") {
    for(i in 1:length(pval)) {
      x <- results[i,2:3]
      if(sum(x) > 0)
        pval[i] <- chisq.test(x,p=c(0.5,0.5))$p.value
    }
    results <- cbind(results, P.value=pval)
  }
  else if(type=="f2") {
    # determine whether marker is autosomal or X-linked
    mar.type <- unlist(sapply(cross$geno,function(a) rep(class(a),ncol(a$data))))

    for(i in 1:length(pval)) {
      if(mar.type[i] == "A") {
        x <- results[i,2:4]
        y <- results[i,5:6]
        if(sum(x) > 0 && sum(y)==0)
          pval[i] <- chisq.test(x,p=c(0.25,0.5,0.25))$p.value
      }
      else {
        x <- results[i,2:3]
        y <- results[i,4:6]
        if(sum(x) > 0 && sum(y)==0)
          pval[i] <- chisq.test(x,p=c(0.5,0.5))$p.value
      }
    }
    results <- cbind(results, P.value=pval)
  }    
  else if(type == "4way") {
    # determine whether marker is autosomal or X-linked
    mar.type <- unlist(sapply(cross$geno,function(a) rep(class(a),ncol(a$data))))

    for(i in 1:length(pval)) {
      if(mar.type[i] == "A") {
        x <- results[i,2:5]
        y <- results[i,-(1:5)]
        if(sum(x) > 0 && sum(y)==0)
          pval[i] <- chisq.test(x,p=c(0.25,0.25,0.25,0.25))$p.value
      }
    }
    results <- cbind(results, P.value=pval)
  }   

  data.frame(chr=rep(names(cross$geno),nmar(cross)),results)
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

#  for(i in 1:nchr(cross)) 
#    cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,
#                                 cross$geno[[i]]$map)
  cross <- makeSSmap(cross)

  cross
}

######################################################################
#
# switch.order: change the marker order on a given chromosome to some
#               specified order
#
######################################################################

switch.order <-
function(cross, chr, order, error.prob=0,
         map.function=c("haldane","kosambi","c-f","morgan"))
{
  map.function <- match.arg(map.function)
  
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

  # save recombination fractions
  flag <- 0
  if(!is.na(match("rf",names(cross)))) {
    rf <- cross$rf
    # determine column within rec fracs
    oldcols <- cumsum(c(0,n.mar))[chr]+seq(along=order)
    newcols <- cumsum(c(0,n.mar))[chr]+order
    rf[oldcols,] <- rf[newcols,]
    rf[,oldcols] <- rf[,newcols]
    colnames(rf)[oldcols] <- colnames(rf)[newcols]
    rownames(rf)[oldcols] <- rownames(rf)[newcols]
    flag <- 1
  }

  # remove any intermediate calculations (except rec fracs),
  #   as they will no longer be meaningful
  cross <- clean(cross)

  if(!is.matrix(cross$geno[[chr]]$map))
    first <- min(cross$geno[[chr]]$map)
  else 
    first <- apply(cross$geno[[chr]]$map,1,min)

  # re-order markers
  cross$geno[[chr]]$data <- cross$geno[[chr]]$data[,order,drop=FALSE]
  m <- seq(0,by=5,length=ncol(cross$geno[[chr]]$data))
  names(m) <- colnames(cross$geno[[chr]]$data)
  if(is.matrix(cross$geno[[chr]]$map)) 
    cross$geno[[chr]]$map <- rbind(m,m)
  else
    cross$geno[[chr]]$map <- m

  # re-estimate rec fracs for re-ordered chromosome
  if(flag==1) {
    temp <- est.rf(subset(cross, chr=chr))$rf
    rf[oldcols,oldcols] <- temp
    cross$rf <- rf
  }

  # re-estimate map
  newmap <- est.map(subset(cross,chr=chr),
                    error.prob=error.prob, map.function=map.function)

  if(!is.matrix(newmap[[1]]))
     cross$geno[[chr]]$map <- newmap[[1]] + first
  else {
     cross$geno[[chr]]$map[1,] <- newmap[[1]][1,] + first[1]
     cross$geno[[chr]]$map[2,] <- newmap[[1]][2,] + first[2]
     rownames(cross$geno[[chr]]$map) <- NULL
   }

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
        stop("If logical, chr argument must have length ", n.chr)
      chr <- sort((1:n.chr)[chr])
    }
        
    if(is.numeric(chr)) {
      # if all negative numbers, convert to positive
      if(all(chr < 1)) chr <- sort((1:n.chr)[chr])
      else chr <- sort(chr)
        
      if(any(chr < 1 | chr > n.chr))
        stop("Chromosome numbers out of range.")
    }
    else {
      if(any(is.na(match(chr,names(x$geno)))))
        stop("Not all chromosome names found.")
      # convert to numeric
      chr <- sort(match(chr,names(x$geno)))
    }

    if(length(chr) != length(unique(chr))) {
      chr <- unique(chr)
      warning("Dropping duplicate chromosomes")
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
      ind[is.na(ind)] <- FALSE
      if(length(ind) != n.ind) 
        stop("If logical, ind argument must have length ", n.ind)
      ind <- (1:n.ind)[ind]
    }
        
    if(is.numeric(ind)) {
      ind <- ind[!is.na(ind)]

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

#pull.chr <-
#function(cross, chr) {
#  warning("pull.chr is deprecated; use subset.cross.")
#  subset.cross(cross, chr)
#}


######################################################################
#
# c.cross: Combine crosses
#
######################################################################

c.cross <-
function(...)
{
  args <- list(...)

  # if only one cross, just return it
  if(length(args)==1) return(args[[1]])

  if(any(sapply(args, function(a) class(a)[2]) != "cross"))
    stop("All arguments must be cross objects.")

  if(length(unique(sapply(args, nchr))) > 1) 
    stop("All arguments must have the same number of chromosomes.")

  x <- args[[1]]
  chr <- names(x$geno)
  n.mar <- nmar(x)
  marnam <- unlist(lapply(x$geno,function(b) colnames(b$data)))

  for(i in 2:length(args)) {
    y <- args[[i]]
    y.marnam <- unlist(lapply(y$geno, function(b) colnames(b$data)))
    if(chr != names(y$geno) || any(n.mar != nmar(y)) || any(marnam != y.marnam))
      stop("All arguments must have the same chromosomes and markers.")
  }
    
  # get all phenotype names
  phenam <- names(x$pheno)
  for(i in 2:length(args))
    phenam <- c(phenam, names(args[[i]]$pheno))
  phenam <- unique(phenam)

  # form big phenotype matrix
  n.ind <- sapply(args,nind)
  pheno <- matrix(nrow=sum(n.ind),ncol=length(phenam))
  colnames(pheno) <- phenam
  pheno <- as.data.frame(pheno)

  for(i in 1:length(phenam)) {
    phe <- vector("list",length(args))
    for(j in 1:length(args)) {
      o <- match(phenam[i],names(args[[j]]$pheno))
      if(is.na(o)) phe[[j]] <- rep(NA,n.ind[j])
      else phe[[j]] <- args[[j]]$pheno[,o]
    }
    pheno[,i] <- unlist(phe)
  }

  # indicator of which cross
  whichcross <- matrix(0,ncol=length(args),nrow=sum(n.ind))
  colnames(whichcross) <- paste("cross",1:length(args),sep="")
  prev <- 0
  for(i in 1:length(args)) {
    wh <- prev + 1:n.ind[i]
    prev <- prev + n.ind[i]
    whichcross[wh,i] <- 1
  }
  pheno <- cbind(pheno,whichcross)

  # crosses must be all the same, or must be combination of F2 and BC
  classes <- sapply(args,function(a) class(a)[1])
  if(length(unique(classes))==1) {
    allsame <- TRUE
    type <- classes[1]
  }
  else {
    if(any(classes != "bc" & classes != "f2")) 
      stop("Experiments must be either the same type or be bc/f2.")
    allsame <- FALSE
    type <- "f2"
    which <- rep(c(0,1)[match(classes,c("bc","f2"))],n.ind)
    pheno <- cbind(pheno,which)
  }

  x$pheno <- pheno

  # create genotype information
  geno <- x$geno
  for(j in 1:nchr(x)) { # drop extraneous stuff
    geno[[j]] <- list(data=geno[[j]]$data, map=geno[[j]]$map)
    class(geno[[j]]) <- class(x$geno[[j]])
  }
  for(i in 2:length(args)) 
    for(j in 1:nchr(x))
      geno[[j]]$data <- rbind(geno[[j]]$data,args[[i]]$geno[[j]]$data)
  
  # if probs exist in each and all have the same
  #     set up values, keep them
  wh <- sapply(args, function(a) match("prob",names(a$geno[[1]])))
  step <- sapply(args,function(a) attr(a$geno[[1]]$prob,"step"))
  error.prob <- sapply(args,function(a) attr(a$geno[[1]]$prob,"error.prob"))
  off.end <- sapply(args,function(a) attr(a$geno[[1]]$prob,"off.end"))
  map.function <- sapply(args,function(a) attr(a$geno[[1]]$prob,"map.function"))
  if(!any(is.na(wh)) && length(unique(step))==1 &&
     length(unique(error.prob))==1 && length(unique(off.end))==1 &&
     length(unique(map.function))==1) {
    if(allsame) { # all same cross type
      for(j in 1:nchr(x)) {
        geno[[j]]$prob <- array(dim=c(sum(n.ind),dim(x$geno[[j]]$prob)[-1]))
        dimnames(geno[[j]]$prob) <- dimnames(x$geno[[j]]$prob)
        prev <- 0
        for(i in 1:length(args)) {
          wh <- prev + 1:n.ind[i]
          prev <- prev + n.ind[i]
          geno[[j]]$prob[wh,,] <- args[[i]]$geno[[j]]$prob
        }
      }
    }
    else { # mixed F2 and BC
      for(j in 1:nchr(x)) {
        wh <- match("f2",classes)
        geno[[j]]$prob <- array(0,dim=c(sum(n.ind),dim(args[[wh]]$geno[[j]]$prob)[-1]))
        dimnames(geno[[j]]$prob) <- dimnames(args[[wh]]$geno[[j]]$prob)
        prev <- 0
        for(i in 1:length(args)) {
          wh <- prev + 1:n.ind[i]
          prev <- prev + n.ind[i]
          if(classes[i]=="f2") 
            geno[[j]]$prob[wh,,] <- args[[i]]$geno[[j]]$prob
          else # backcross
            geno[[j]]$prob[wh,,1:2] <- args[[i]]$geno[[j]]$prob
        }
      }
    }    

    for(j in 1:nchr(x)) {
      attr(geno[[j]]$prob,"step") <- step[1]
      attr(geno[[j]]$prob,"error.prob") <- error.prob[1]
      attr(geno[[j]]$prob,"off.end") <- off.end[1]
      attr(geno[[j]]$prob,"map.function") <- map.function[1]
    }
  }

  # if draws exist in each and all have the same
  #     set up values, keep them
  wh <- sapply(args, function(a) match("draws",names(a$geno[[1]])))
  step <- sapply(args,function(a) attr(a$geno[[1]]$draws,"step"))
  error.prob <- sapply(args,function(a) attr(a$geno[[1]]$draws,"error.prob"))
  off.end <- sapply(args,function(a) attr(a$geno[[1]]$draws,"off.end"))
  map.function <- sapply(args,function(a) attr(a$geno[[1]]$draws,"map.function"))
  ndraws <- sapply(args,function(a) dim(a$geno[[1]]$draws)[3])
  if(!any(is.na(wh)) && length(unique(step))==1 &&
     length(unique(error.prob))==1 && length(unique(off.end))==1 &&
     length(unique(map.function))==1 && length(unique(ndraws))==1) {
    for(j in 1:nchr(x)) {
      geno[[j]]$draws <- array(0,dim=c(sum(n.ind),dim(x$geno[[j]]$draws)[-1]))
      dimnames(geno[[j]]$draws) <- dimnames(x$geno[[j]]$draws)
      prev <- 0
      for(i in 1:length(args)) {
        wh <- prev + 1:n.ind[i]
        prev <- prev + n.ind[i]
        geno[[j]]$draws[wh,,] <- args[[i]]$geno[[j]]$draws
      }

      attr(geno[[j]]$draws,"step") <- step[1]
      attr(geno[[j]]$draws,"error.prob") <- error.prob[1]
      attr(geno[[j]]$draws,"off.end") <- off.end[1]
      attr(geno[[j]]$draws,"map.function") <- map.function[1]
    }
  }

  x <- list(geno=geno, pheno=pheno)
  class(x) <- c(type,"cross")
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
  chrtype <- sapply(cross$geno, class)

  # drop individuals whose sex or pgm is missing if X chr is included
  if(any(chrtype=="X")) {
    sexpgm <- getsex(cross)
    keep <- rep(TRUE,nind(cross))
    flag <- 0
    if(!is.null(sexpgm$sex)) {
      if(any(is.na(sexpgm$sex))) {
        keep[is.na(sexpgm$sex)] <- FALSE
        flag <- 1
      }
    }
    if(!is.null(sexpgm$pgm)) {
      if(any(is.na(sexpgm$pgm))) {
        keep[is.na(sexpgm$pgm)] <- FALSE
        flag <- 1
      }
    }
    if(flag) {
      warning("Dropping ", sum(!keep), " individuals with missing sex or pgm.\n")
      cross <- subset(cross, ind=keep)
      if(!is.null(addcovar)) {
        if(!is.matrix(addcovar)) addcovar <- addcovar[keep]
        else addcovar <- addcovar[keep,]
      }
      if(!is.null(intcovar)) {
        if(!is.matrix(intcovar)) intcovar <- intcovar[keep]
        else intcovar <- intcovar[keep,]
      }
    }
  }

  # check phenotypes
  if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
  if(pheno.col < 1 || pheno.col > nphe(cross))
    stop("Specified phenotype column is invalid.")
  if(!is.numeric(cross$pheno[,pheno.col]))
    stop("Chosen phenotype is not numeric, and needs to be.")

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
      warning("Dropping ", sum(wh), " individuals with missing covariates.\n")
    }
  }

  # make sure columns of intcovar are contained in addcovar
  if(!is.null(intcovar)) {
    if(is.null(addcovar)) {
      addcovar <- intcovar
      n.addcovar <- n.intcovar
      warning("addcovar forced to contain all columns of intcovar\n")
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
function(r, type=c("self","sib"), chrtype=c("A","X"), expand=TRUE)
{
  # type of RI lines
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)

  if(type=="self") {
    if(expand) return(r*2/(1+2*r))
    else return(r/2/(1-r))
  }
  else {
    if(chrtype=="A") { # autosome / sib mating
      if(expand) return(r*4/(1+6*r))
      else return(r/(4-6*r))
    }
    else { # X chromosome/ sib mating
      if(expand) return(8/3*r/(1+4*r))
      else return(3/8*r/(1-1.5*r))
    }
  }
}

######################################################################
# pull.geno
######################################################################
pull.geno <-
function(cross)
{
  X <- cross$geno[[1]]$data
  if(nchr(cross) > 1)
    for(i in 2:nchr(cross))
      X <- cbind(X, cross$geno[[i]]$data)
  X
}

######################################################################
# lodint: function to get lod support interval
######################################################################
lodint <-
function(results, chr, drop=1.5)
{
  results <- results[results[,1]==chr,]

  if(all(is.na(results[,3]))) return(NULL)

  maxlod <- max(results[,3],na.rm=TRUE)
  w <- which(!is.na(results[,3]) & results[,3] == maxlod)
  o <- range(which(!is.na(results[,3]) & results[,3] > maxlod-drop))

  if(length(o)==0) o <- c(1,nrow(results))

  else {
    if(o[1] > 1) o[1] <- o[1]-1
    if(o[2] < nrow(results)) o[2] <- o[2]+1
  }

  results <- results[c(o[1],w,o[2]),]
  class(results) <- c("scanone","data.frame")

  results
}


######################################################################
# bayesint: function to get Bayesian probability interval
######################################################################
bayesint <-
function(results, chr, prob=0.95)
{
  results <- results[results[,1]==chr,]

  if(all(is.na(results[,3]))) return(NULL)

  loc <- results[,2]
  width <- diff(( c(loc[1],loc) + c(loc, loc[length(loc)]) )/ 2)

  area <- 10^results[,3]*width
  area <- area/sum(area)

  o <- rev(order(results[,3]))

  cs <- o
  for(i in 1:length(o))
    cs[i] <- sum(area[min(o[1:i]):max(o[1:i])])

  wh <- min((1:length(loc))[cs >= prob])
  int <- range(o[1:wh])

  results[c(int[1],o[1],int[2]),]
}


  
######################################################################
# makeSSmap: convert a genetic map, or the genetic maps in a cross
#            object, to be sex-specific (i.e., 2-row matrices)
######################################################################
makeSSmap <-
function(cross)
{
  if(class(cross)[1] == "map") {
    # input object is a genetic map
    for(i in 1:length(cross)) {
      if(!is.matrix(cross[[i]]))
        cross[[i]] <- rbind(cross[[i]], cross[[i]])
    }
  }
  else { # input object is assumed to be a "cross" object
    n.chr <- nchr(cross)
    for(i in 1:n.chr) {
      if(!is.matrix(cross$geno[[i]]$map))
        cross$geno[[i]]$map <-
          rbind(cross$geno[[i]]$map, cross$geno[[i]]$map)
    }
  }

  cross
}

######################################################################
# comparecrosses: verify that two cross objects have identical
#                 classes, chromosomes, markers, genotypes, maps,
#                 and phenotypes
######################################################################
comparecrosses <-
function(cross1, cross2, tol=1e-5)
{
  # both are of class "cross"
  if(is.na(match("cross", class(cross1)))) 
    stop("cross1 is not a cross object.")
  if(is.na(match("cross", class(cross2)))) 
    stop("cross2 is not a cross object.")

  # classes are the same
  if(any(class(cross1) != class(cross2))) 
    stop("crosses are not the same type.")

  if(nchr(cross1) != nchr(cross2)) 
    stop("crosses do not have the same number of chromosomes.")

  if(any(names(cross1$geno) != names(cross2$geno)))
    stop("Chromosome names do not match.")

  if(any(nmar(cross1) != nmar(cross2)))
    stop("Number of markers per chromosome do not match.")

  mnames1 <- unlist(lapply(cross1$geno, function(a) colnames(a$data)))
  mnames2 <- unlist(lapply(cross2$geno, function(a) colnames(a$data)))
  if(any(mnames1 != mnames2)) {
#    stop("Markers names do not match.")
    for(i in 1:nchr(cross1)) 
      if(any(colnames(cross1$geno[[i]]$data) != colnames(cross2$geno[[i]]$data)))
        stop("Marker names on chr ", names(cross1$geno)[i], " don't match.")
  }


  chrtype1 <- sapply(cross1$geno, class)
  chrtype2 <- sapply(cross2$geno, class)
  if(any(chrtype1 != chrtype2))
    stop("Chromosome types (autosomal vs X) do not match.")

  for(i in 1:nchr(cross1)) {
    if(any(abs(diff(cross1$geno[[i]]$map) - diff(cross2$geno[[i]]$map)) > tol))
      stop("Genetic maps for chromosome ", names(cross1$geno)[i],
           " do not match.")
      
    if(abs(cross1$geno[[i]]$map[1] - cross2$geno[[i]]$map[1]) > tol)
      warning("Initial marker positions for chromosome ", names(cross1$geno)[i],
              " do not match.")
  }

  if(nind(cross1) != nind(cross2))
    stop("Number of individuals do not match.")

  for(i in 1:nchr(cross1)) {
    g1 <- cross1$geno[[i]]$data
    g2 <- cross2$geno[[i]]$data
    if(any((is.na(g1) & !is.na(g2)) | (!is.na(g1) & is.na(g2)) |
           (!is.na(g1) & !is.na(g2) & g1!=g2))) 
      stop("Genotype data for chromosome ", names(cross1$geno)[i],
           " do not match.")
  }

  if(nphe(cross1) != nphe(cross2))
    stop("Number of phenotypes do not match.")

  if(any(names(cross1$pheno) != names(cross2$pheno)))
    stop("Phenotype names do not match.")

  for(i in 1:nphe(cross1)) {
    phe1 <- cross1$pheno[,i]
    phe2 <- cross2$pheno[,i]
    if(is.numeric(phe1) & is.numeric(phe2)) {
      if(any((is.na(phe1) & !is.na(phe2)) | (!is.na(phe1) & is.na(phe2)) |
             (!is.na(phe1) & !is.na(phe2) & abs(phe1-phe2) > tol))) {
        stop("Data for phenotype ", names(cross1$pheno)[i],
             " do not match.")
      }
    }
    else {
      if(any((is.na(phe1) & !is.na(phe2)) | (!is.na(phe1) & is.na(phe2)) |
             (!is.na(phe1) & !is.na(phe2) &
              as.character(phe1) != as.character(phe2)))) {
        stop("Data for phenotype ", names(cross1$pheno)[i], " do not match.")
      }
    }
  }

  cat("\tCrosses are identical.\n")
}


######################################################################
# move marker
# Move a marker to a new chromosome...placed at the end
######################################################################
movemarker <-
function(cross, marker, newchr, newpos)
{
  mnames <- unlist(lapply(cross$geno,function(a) colnames(a$data)))
  chr <- rep(names(cross$geno),nmar(cross))
  pos <- unlist(lapply(cross$geno,function(a) 1:ncol(a$data)))
  oldindex <- match(marker, mnames)

  # Marker found precisely once?
  if(is.na(oldindex)) stop(marker, " not found.\n")
  if(length(oldindex) > 1) stop(marker, " found multiple times.\n")

  if(is.na(match(newchr,names(cross$geno))))
    stop("Chromosome ", newchr, " not found.\n")

  chr <- chr[oldindex]
  pos <- pos[oldindex]

  # pull out genotype data
  g <- cross$geno[[chr]]$data[,pos]

  # delete marker
  if(nmar(cross)[chr] == 1)  { # only marker on that chromosome, so drop the chromosome
    cross$geno <- cross$geno[-match(chr,names(cross$geno))]
    delchr <- TRUE
  }
  else {
    delchr <- FALSE
    cross$geno[[chr]]$data <- cross$geno[[chr]]$data[,-pos,drop=FALSE]
    if(is.matrix(cross$geno[[chr]]$map))
      cross$geno[[chr]]$map <- cross$geno[[chr]]$map[,-pos,drop=FALSE]
    else
      cross$geno[[chr]]$map <- cross$geno[[chr]]$map[-pos]
  }

  if(missing(newpos)) {
    # add marker to end of new chromosome
    n.mar <- nmar(cross)[newchr]
    cross$geno[[newchr]]$data <- cbind(cross$geno[[newchr]]$data,g)
    colnames(cross$geno[[newchr]]$data)[n.mar+1] <- marker
  
    if(is.matrix(cross$geno[[newchr]]$map)) {
      cross$geno[[newchr]]$map <- cbind(cross$geno[[newchr]]$map,
                                        cross$geno[[newchr]]$map[,n.mar]+10)
      colnames(cross$geno[[newchr]]$map)[n.mar+1] <- marker
    }
    else {
      cross$geno[[newchr]]$map <- c(cross$geno[[newchr]]$map,
                                    cross$geno[[newchr]]$map[n.mar]+10)
      names(cross$geno[[newchr]]$map)[n.mar+1] <- marker
    }
  }
  else {
    # add marker to the specified position
    dat <- cross$geno[[newchr]]$data
    map <- cross$geno[[newchr]]$map

    if(length(newpos) != 1)
      stop("newpos should be a single number.")

    if(is.matrix(map)) { # sex-specific maps
      wh <- which(map[1,] < newpos)
      if(length(wh) == 0) { # place in first spot
        map <- cbind(c(newpos,map[2,1]-(map[1,1]-newpos)),map)
        colnames(map)[1] <- marker
      }
      else {
        wh <- max(wh)
        if(wh == ncol(map)) { # place at end of chromosome
          map <- cbind(map,c(newpos,map[2,ncol(map)]+(newpos-map[1,ncol(map)])))
          colnames(map)[ncol(map)] <- marker
        }
        else {
          left <- map[,wh]
          right <- map[,wh+1]
          newpos2 <- (newpos-left[1])/(right[1]-left[1])*(right[2]-left[2])+left[2]
          map <- cbind(map[,1:wh], c(newpos,newpos2), map[,-(1:wh)])
          colnames(map)[wh+1] <- marker
        }
      }
    }
    else {
      wh <- which(map < newpos)
      if(length(wh) == 0) { # place in first position
        map <- c(newpos,map)
        names(map)[1] <- marker
      }
      else {
        wh <- max(wh)
        if(wh == length(map)) { # place in last position
          map <- c(map,newpos)
          names(map)[length(map)] <- marker
        }
        else {
          map <- c(map[1:wh],newpos,map[-(1:wh)])
          names(map)[wh+1] <- marker
        }
      }
    }
    cross$geno[[newchr]]$map <- map

    if(length(wh)==0) { # place marker in first position
      dat <- cbind(g, dat)
      colnames(dat)[1] <- marker
    }
    else if(wh == ncol(dat)) { # place marker in last position
      dat <- cbind(dat, g)
      colnames(dat)[ncol(dat)] <- marker
    }
    else { # place marker in the middle
      dat <- cbind(dat[,1:wh],g,dat[,-(1:wh)])
      colnames(dat)[wh+1] <- marker
    }
    cross$geno[[newchr]]$data <- dat
    
    # make sure the marker names for the data and the genetic map match
    colnames(cross$geno[[newchr]]$data) <- names(cross$geno[[newchr]]$map)
  }

  # update genoprob, errorlod, argmax, draws, rf
  if(!is.na(match("rf",names(cross)))) {
    # reorder the recombination fractions
    #  -- a bit of pain, 'cause we need LODs in upper triangle
    #     and rec fracs in lower triangle
    newmar <- unlist(lapply(cross$geno,function(a) colnames(a$data)))

    rf <- cross$rf
    lods <- rf;lods[lower.tri(rf)] <- t(rf)[lower.tri(rf)]
    rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]
    lods <- lods[newmar,newmar]
    rf <- rf[newmar,newmar]
    rf[upper.tri(rf)] <- lods[upper.tri(rf)]

    cross$rf <- rf
  }

  if(!delchr) thechr <- c(chr,newchr)
  else thechr <- newchr
  for(i in thechr) {
    tempg <- cross$geno[[i]]
    tempx <- subset(cross, chr=i)
    if(!is.na(match("prob", names(tempg))))
      atp <- attributes(tempg$prob) 
    if(!is.na(match("draws", names(tempg)))) {
      at <- attributes(listeria$geno[[1]]$draws)
      tempg$draws <- sim.geno(tempx,
                              n.draws=at$dim[3],
                              step=at$step,
                              off.end=at$off.end,
                              map.function=at$map.function,
                              error.prob=at$error.prob)$geno[[1]]$draws
    }
    if(!is.na(match("argmax",names(tempg)))) {
      at <- attributes(listeria$geno[[1]]$argmax)
      tempg$argmax <- argmax.geno(tempx,
                                 step=at$step,
                                 off.end=at$off.end,
                                 map.function=at$map.function,
                                 error.prob=at$error.prob)$geno[[1]]$argmax
    }
    if(!is.na(match("errorlod",names(tempg)))) {
      at <- attributes(listeria$geno[[1]]$errorlod)
      tempg$errorlod <- argmax.geno(tempx,
                                    map.function=at$map.function,
                                    error.prob=at$error.prob)$geno[[1]]$errorlod
    }
    if(!is.na(match("prob",names(tempg)))) 
      tempg$prob <- calc.genoprob(tempx,
                                  step=atp$step,
                                  off.end=atp$off.end,
                                  map.function=atp$map.function,
                                  error.prob=atp$error.prob)$geno[[1]]$prob

    cross$geno[[i]] <- tempg
  }

  cross
}

######################################################################
#
# summary.map
#
# Give a short summary of a genetic map object.
# 
######################################################################
summary.map <- 
function(object, ...)
{
  map <- object
  if(length(class(map))>1 && class(map)[2] == "cross") # a cross object
    map <- pull.map(map)
  
  n.chr <- length(map)
  chrnames <- names(map)
  if(is.matrix(map[[1]])) { # sex-specific map
    sexsp <- TRUE
    n.mar <- sapply(map,ncol)
    tot.mar <- sum(n.mar)
    fmap <- lapply(map,function(a) a[1,])
    mmap <- lapply(map,function(a) a[2,])

    len.f <- sapply(fmap,function(a) diff(range(a)))
    len.m <- sapply(mmap,function(a) diff(range(a)))
    avesp.f <- sapply(fmap,function(a) mean(diff(a)))
    avesp.m <- sapply(mmap,function(a) mean(diff(a)))
    totlen.f <- sum(len.f)
    totlen.m <- sum(len.m)

    tot.avesp.f <- mean(unlist(lapply(fmap,diff)))
    tot.avesp.m <- mean(unlist(lapply(mmap,diff)))
                    
    output <- rbind(cbind(n.mar,len.f,len.m,avesp.f,avesp.m),
                    c(tot.mar,totlen.f,totlen.m,tot.avesp.f,tot.avesp.m))
    dimnames(output) <- list(c(chrnames,"overall"),
                             c("n.mar","length.female","length.male",
                               "ave.spacing.female","ave.spacing.male"))
  }                   
  else {
    sexsp=FALSE
    n.mar <- sapply(map,length)
    len <- sapply(map,function(a) diff(range(a)))
    tot.mar <- sum(n.mar)

    len <- sapply(map,function(a) diff(range(a)))
    avesp <- sapply(map,function(a) mean(diff(a)))
    totlen <- sum(len)
    tot.avesp <- mean(unlist(lapply(map,diff)))
                    
    output <- rbind(cbind(n.mar,len,avesp),
                    c(tot.mar,totlen,tot.avesp))
    dimnames(output) <- list(c(chrnames,"overall"),
                             c("n.mar","length","ave.spacing"))

  }

  output <- list(summarytable=output,sexsp=sexsp)
  class(output) <- "summary.map"
  output
}


######################################################################
#
# print.summary.map
#
# Print out the result of summary.map()
# 
######################################################################
print.summary.map <-
function(x, ...)  
{
  if(x[[2]]) cat("Sex-specific map\n\n")
  else cat("Sex-averaged map\n\n")

  x <- x[[1]]
  x <- apply(x,2,round,1)
  print(x)
}
  
######################################################################
#
# convert.scanone
#
# Convert scanone output from the format for R/qtl ver 0.97 to
# that for R/qtl ver 0.98
# (previously, inter-maker locations named loc*.c*; now c*.loc*)
#
######################################################################
convert.scanone <-
function(output)
{  
  rn <- rownames(output)
  o <- grep("^loc\-*[0-9]+(\.[0-9]+)*\.c[0-9A-Za-z]+$", rn)
  if(length(o) > 0) {
    temp <- rn[o]
    temp <- strsplit(temp,"\\.")
    temp <- sapply(temp, function(a)
                   paste(a[c(length(a),1:(length(a)-1))],collapse="."))
    rownames(output)[o] <- temp
  }
  output
}
            
######################################################################
# find.pheno
#
# utility to get pheno number given pheno name
######################################################################
find.pheno <-
function( cross,  pheno )
  seq( ncol( cross$pheno ))[match(pheno,names(cross$pheno))]

######################################################################
# find.flanking
#
# utility to get flanking and/or closest marker to chr and pos
######################################################################
find.flanking <-
function( cross, chr, pos)
{
  map = pull.map(cross)

  if(is.matrix(map[[1]]) && nrow(map[[1]]) > 1) 
    stop("This function works only for crosses with sex-averaged maps.")

  if(length(chr) == 1 && length(pos) > 1) {
    chr <- rep(chr,length(pos))
  }

  marker = NULL
  for (i in seq(length(chr))) {
    tmp = map[[chr[i]]]-pos[i]
    m = names(map[[chr[i]]])
    left = sum(tmp < 0)
    at = sum(tmp == 0)
    right = sum(tmp > 0)
    f <- if (at > 0)
      left+at[c(1,length(at))]
    else {
      if (right > 0)
        c(left,left+at+1)
      else
        c(left,left+at)
    }
    marker = rbind(marker,m[f[c(1:2,order(abs(tmp[f]))[1])]])
  }
  dimnames(marker) <- list(paste("chr",chr,":",pos,sep=""),
                           c("left","right","close"))
  as.data.frame(marker)
}

######################################################################
# strip.partials
#
# Removes partially informative genotypes in an intercross.
#
# Input: a cross object; if verbose=TRUE, a statement regarding the
#        number of genotypes removed is printed.
######################################################################
strip.partials <-
function(cross, verbose=TRUE)
{
  type <- class(cross)[1]
  if(type != "f2") 
    stop("This is for intercrosses only")

  n.removed <- 0
  for(i in 1:nchr(cross)) {
    g <- cross$geno[[i]]$data
    wh <- !is.na(g) & g>3
    if(any(wh)) {
      g[wh] <- NA
      cross$geno[[i]]$data <- g
      n.removed <- n.removed + sum(wh)
    }
  }
  if(verbose) {
    if(n.removed == 0) cat(" --Didn't remove any genotypes.\n")
    else cat("Removed", n.removed, "genotypes.\n")
  }
  cross
}

######################################################################
# comparegeno
######################################################################
comparegeno <-
function(cross, what=c("proportion","number"))
{
  what <- match.arg(what)
  g <- pull.geno(cross)
  g[is.na(g)] <- 0
  n.ind <- nrow(g)
  n.mar <- ncol(g)
  z <- .C("R_comparegeno",
          as.integer(g),
          as.integer(n.ind),
          as.integer(n.mar),
          n.match=as.integer(rep(0,n.ind^2)),
          n.missing=as.integer(rep(0,n.ind^2)),
          PACKAGE="qtl")

  if(what=="number") return(matrix(z$n.match,n.ind,n.ind))
  else return(matrix(z$n.match/(n.mar-z$n.missing),n.ind,n.ind))
}
  
          


# end of util.R


