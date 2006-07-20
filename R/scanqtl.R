######################################################################
#
# scanqtl.R
#
# copyright (c) 2002-6, Hao Wu, The Jackson Laboratory
#                       and Karl W. Broman, Johns Hopkins University
# last modified Jul, 2006
# first written Apr, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: scanqtl
#
######################################################################

scanqtl <-
  function(cross, pheno.col=1, chr, pos, covar=NULL, formula, method=c("imp"),
           incl.markers=FALSE, verbose=TRUE)
{
  type <- class(cross)[1]
  chrtype <- sapply(cross$geno,class)
  sexpgm <- getsex(cross)
  
  # input data checking
  if( !sum(class(cross) == "cross") )
    stop("The first input variable must be  an object of class cross")
  if( length(chr) != length(pos))
    stop("Input chr and pos must have the same length")
  # note that input chr is a vector and pos is a list

  # check the input covariate, if any
  if(!is.null(covar)) {
    if(nrow(covar) != nind(cross))
      stop("Input covariate has wrong size")
  }
  # check the input pheno.col
  if(pheno.col>ncol(cross$pheno))
    stop("Wrong phenotype column number")
  
  method <- match.arg(method)

  ichr <- match(chr, names(cross$geno))
  if(any(is.na(ichr))) {
    err <- paste("There's no chromosome number ", chr[is.na(ichr)],
                 "in input cross object")
    stop(err)
  }

  # if formula is missing, make one.
  # All QTLs and covariates will be additive by default
  n.qtl <- length(chr)
  n.covar <- length(covar)
  if(missing(formula)) {
    tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
    formula <- "y~Q1"
    if(n.qtl > 1) 
      for (i in 2:n.qtl) 
        formula <- paste(formula, tmp.Q[i], sep="+")
    if (n.covar) { # if covariate is not empty
      tmp.C <- names(covar) # covariate term names
      for(i in 1:n.covar)
        formula <- paste(formula, tmp.C[i], sep="+")
    }
    formula <- as.formula(formula)
  }
  else {
    # include all input QTLs and covariates in the formula additively
    formula.str <- deparse(formula) # deparse formula as a string
    for(i in 1:n.qtl) { # loop thru the QTLs
      qtl.term <- paste("Q", i, sep="")
      if( length(grep(qtl.term, formula.str, ignore.case=TRUE))==0 )
        # this term is not in the formula
        # add it to the formula
        formula.str <- paste(formula.str, qtl.term, sep="+")
    }
    if(n.covar) { # covariate is not empty
      for(i in 1:n.covar) {
        covar.term <- names(covar)[i]
        if( length(grep(covar.term, formula.str, ignore.case=TRUE))==0 )
        # this term is not in the formula
        # add it to the formula
          formula.str <- paste(formula.str, covar.term, sep="+")
      }
    }
    formula <- as.formula(formula.str)
  }
  
  # find the chromosome with multiple QTLs
  # indices for chromosomes with multiple QTLs
  idx.varied <- NULL
  indices <- pos  ## added by Karl 8/23/05
  for(i in 1:length(pos)) {
    l <- length(pos[[i]] )
    if( l >= 2 ) {
      # if there're more than two elements in pos, issue warning message
      if(l > 2) {
        msg <- "There are more than two elements in "
        msg <- paste(msg, i, "th input pos.")
        msg <- paste(msg, "The first two are taken as starting and ending position.")
        warning(msg)
      }

      # user specified a range
      # find all markers in this range
      idx.varied <- c(idx.varied, i) 
      # make the genetic map on this chromosome
      if(!("draws" %in% names(cross$geno[[1]]))) # there's no draw in input cross object
        stop("You need to first run sim.geno().")
      # make genetic map
      if("map" %in% names(attributes(cross$geno[[ichr[i]]]$draws)))
        map <- attr(cross$geno[[ichr[i]]]$draws,"map")
      else {
        stp <- attr(cross$geno[[ichr[i]]]$draws, "step")
        oe <- attr(cross$geno[[ichr[i]]]$draws, "off.end")
      
        if("stepwidth" %in% names(attributes(cross$geno[[ichr[i]]]$draws)))
          stpw <- attr(cross$geno[[ichr[i]]]$draws, "stepwidth")
        else stpw <- "fixed"
        map <- create.map(cross$geno[[ichr[i]]]$map,stp,oe,stpw)
      }
      # pull out the female map if there are sex-specific maps
      if(is.matrix(map)) map <- map[1,]

      indices[[i]] <- seq(along=map)
      if(!incl.markers) { # equally spaced positions
        step <- attr(cross$geno[[ichr[i]]]$draws,"step")
        eq.sp.pos <- seq(min(map), max(map), by=step)
        wh.eq.pos <- match(eq.sp.pos, map)
        map <- map[wh.eq.pos]
        indices[[i]] <- indices[[i]][wh.eq.pos]
      }

      # locate the markers given starting and ending postion
      # we should do this before or after incl.markers?
      start <- pos[[i]][1]
      end <- pos[[i]][2]
      # replace pos[[i]] (a range) by the marker positions within the range
      # extend the position to the nearest markers outside the ranges
      tmp <- which( (map - start)<0 )
      if(length(tmp) != 0) # starting position is after the first marker
        start <- map[max(tmp)]
      tmp <- which( (end-map) < 0 )
      if(length(tmp) != 0) # ending position is before the last marker
        end <- map[min(tmp)]
      pos[[i]] <- as.vector( map[(map>=start)&(map<=end)] )
      indices[[i]] <- indices[[i]][(map>=start)&(map<=end)]
    }
    else { # fixed position rather than range
# Hao asked me to comment these two lines out      
#      pos[[i]] <- locatemarker(cross$geno[[as.character(chr[i])]]$map,
#                               pos[[i]], chr[i], "draws")
    }
  }
  # Now, pos contains all the marker positions for all chromosomes
                  
  #########################
  # Now start general scan
  #########################
  # There might be several chromosomes with multiple QTLs
  # Use one loop
  
  # number of chromosomes with multiple positions to be scanned
  n.idx.varied <- length(idx.varied) 
  n.loop <- 1 # total number of loops
  if(n.idx.varied != 0) { # there IS some chromosomes with multiple QTL
    # vector to indicate the positions indices for those chromosomes
    idx.pos <- rep(0, n.idx.varied)
    l.varied <- NULL
    for(i in 1:n.idx.varied) {
      l.varied[i] <- length(pos[[idx.varied[i]]])
      n.loop <- n.loop * l.varied[i]
    }
    # initialize output variable
    result <- array(rep(0, n.loop), rev(l.varied))
  }
  else { # fixed QTL model (no scanning)
    qtl <- makeqtl(cross, chr=chr, pos=unlist(pos))
    result <- fitqtl(cross$pheno[,pheno.col], qtl, covar=covar,
                     formula=formula, method=method, dropone=FALSE,
                     get.ests=FALSE)
    result <- result[1]
    names(result) <- "LOD"
    class(result) <- "scanqtl"
    attr(result, "method") <- method
    attr(result, "formula") <- formula
    return(result)
  }

  # loop thru all varied QTLs
  if(verbose) cat(" ",n.loop, "models to fit\n")
  current.pos <- NULL ## added by Karl 8/23/05
  for(i in 1:n.loop) {
    # find the indices for positions
    remain <- i
    if(n.idx.varied > 1) {
      for(j in 1:(n.idx.varied-1)) {
        ns <- 1
        for( k in (j+1):n.idx.varied )
          ns <- ns * length(pos[[idx.varied[k]]])
        idx.pos[j] <- floor(remain / ns) + 1
        remain <- remain - (idx.pos[j]-1) * ns
        # remain cannot be zero
        if(remain == 0) {
          idx.pos[j] <- idx.pos[j] - 1
          remain <- remain + ns
        }
      }
    }
    idx.pos[n.idx.varied] <- remain

    # make an QTL object 
    pos.tmp <- NULL
    for(j in 1:length(pos)) {
      if(j %in% idx.varied) {
        idx.tmp <- which(j==idx.varied)
        pos.tmp <- c(pos.tmp, pos[[j]][idx.pos[idx.tmp]])
      }
      else
        pos.tmp <- c(pos.tmp, pos[[j]])
    }

    # this bit revised by Karl 8/23/05; now we make the qtl object
    #     once, and copy stuff over otherwise
    if(is.null(current.pos)) {
      qtl.obj <- makeqtl(cross, chr, pos.tmp)
      current.pos <- pos.tmp
    }
    else {
      for(kk in seq(along=pos.tmp)) {
        if(pos.tmp[kk] != current.pos[kk]) {
          u <- abs(pos.tmp[kk]-pos[[kk]])
          w <- indices[[kk]][u==min(u)]
          qtl.obj$geno[,kk,] <- cross$geno[[ichr[kk]]]$draws[,w,]

          if(chrtype[ichr[kk]]=="X")
            qtl.obj$geno[,kk,] <-
              reviseXdata(type,"full",sexpgm,draws=qtl.obj$geno[,kk,,drop=FALSE],
                          cross.attr=attributes(cross))

          current.pos <- pos.tmp
        }
      }
    }
    # end of Karl's 8/23/05 addition

    # fit QTL, don't do drop one at a time
    fit <- fitqtl(cross$pheno[,pheno.col], qtl=qtl.obj, covar=covar,
                  formula=formula, method=method, dropone=FALSE,
                  get.ests=FALSE)
  
    if(verbose) { # feedback to let you know what's happening
      if(n.loop < 10 || 
         (n.loop < 40 && i/4==round(i/4)) ||
         (n.loop < 250 && i==round(i,-1)) ||
         i/20==round(i/20))
        cat("    ", i,"/", n.loop, "\n")
    }
      
    # assign to result matrix
    #     Note: [[1]][1,4] picks out the LOD score 
    result[i] <- fit[[1]][1,4]
  }

  # make the row and column names for the result matrix
  dnames <- list(NULL)
  for(i in 1:n.idx.varied) {
    i.chr <- chr[idx.varied[n.idx.varied-i+1]]
    i.pos <- pos[[idx.varied[n.idx.varied-i+1]]]
    dnames[[i]] <- paste( paste("Chr", i.chr,sep=""),
                            i.pos, sep="@")
  }
  dimnames(result) <- dnames
  
  class(result) <- "scanqtl"
  attr(result, "method") <- method
  attr(result, "formula") <- formula
  result
}


#summary.scanqtl <- function(object, ...)
#{
#}

#print.summary.qtl <- function(x, ...)
#{
#}

# end of scanqtl.R
