######################################################################
#
# vbscan.R
#
# copyright (c) 2001-2, Karl W Broman, Johns Hopkins University
# last modified June, 2002
# first written May, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: vbscan
#
######################################################################

######################################################################
#
# vbscan: scan genome for a quantitative phenotype for which some
# individuals' phenotype is undefined (for example, the size of a
# lesion, where some individuals have no lesion).
#
######################################################################

vbscan <-
function(cross, pheno.col=1, x.treatment=c("simple","full"),
         upper=FALSE, method="em", maxit=4000, tol=1e-4)
{
  method <- match.arg(method)
  type <- class(cross)[1]
  x.treatment <- match.arg(x.treatment)

  # check arguments are okay
  if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
  if(pheno.col > nphe(cross))
    stop("Specified phenotype column exceeds the number of phenotypes")
  y <- cross$pheno[,pheno.col]

  # remove individuals with missing phenotypes
  if(any(is.na(y))) {
    drop <- (seq(along=y))[is.na(y)]
    y <- y[-drop]
    for(i in 1:length(cross$geno)) 
      cross$geno[[i]]$prob <- cross$geno[[i]]$prob[-drop,,]
  }

  # modify phenotypes
  if(upper) {
    if(!any(y == Inf)) y[y==max(y)] <- Inf
  }
  else {
    if(!any(y == -Inf)) y[y==min(y)] <- -Inf
  }
  survived <- rep(0,length(y))
  survived[y == -Inf | y == Inf] <- 1

  # The following line is included since .C() doesn't accept Infs
  y[y == -Inf | y == Inf] <- 99999

  n.chr <- nchr(cross)
  results <- NULL

  for(i in 1:n.chr) {
    # make sure inferred genotypes or genotype probabilities are available
    if(is.na(match("prob",names(cross$geno[[i]])))) {
      cat(" -Calculating genotype probabilities\n")
      cross <- calc.genoprob(cross)
    }

    n.pos <- dim(cross$geno[[i]]$prob)[2]
    n.ind <- length(y)

    chrtype <- class(cross$geno[[i]])
    if(chrtype=="X") sexpgm <- getsex(cross)
    else sexpgm <- NULL

    gen.names <- getgenonames(type,chrtype,x.treatment,sexpgm)
    n.gen <- length(gen.names)

    # Update X chromosome
    if(chrtype=="X" && (type=="f2" || type=="f2ss" || type=="bc"))
      cross$geno[[i]]$prob <- fixXdata(type, x.treatment, sexpgm,
                                       prob=cross$geno[[i]]$prob)

    z <- .C("R_vbscan",
            as.integer(n.pos),
            as.integer(n.ind),
            as.integer(n.gen),
            as.double(cross$geno[[i]]$prob),
            as.double(y),
            as.integer(survived),
            lod=as.double(rep(0,(4+2*n.gen)*n.pos)),
            as.integer(maxit),
            as.double(tol))

    map <- create.map(cross$geno[[i]]$map,
                      attr(cross$geno[[i]]$prob,"step"),
                      attr(cross$geno[[i]]$prob,"off.end"))
    if(is.matrix(map)) map <- map[1,]

    res <- data.frame(chr=rep(names(cross$geno)[i],length(map)),
                      pos = map,
                      matrix(z$lod,nrow=n.pos,byrow=TRUE))

    w <- names(map)
    o <- grep("^loc\-*[0-9]+",w)
    if(length(o) > 0)
      w[o] <- paste(w[o],names(cross$geno)[i],sep=".c")
    rownames(res) <- w
    
    colnames(res) <- c("chr","pos","lod","lod.p","lod.mu",
                       paste("pi",gen.names,sep="."),
                       paste("mu",gen.names,sep="."), "sigma")

    z <- res

    # if different number of columns from other chromosomes,
    #     expand to match
    if(!is.null(results) && ncol(z) != ncol(results)) {
      cnz <- colnames(z)
      cnr <- colnames(results)
      wh.zr <- match(cnz,cnr)
      wh.rz <- match(cnr,cnz)
      if(all(!is.na(wh.rz))) {
        newresults <- data.frame(matrix(NA,nrow=nrow(results),ncol=ncol(z)))
        dimnames(newresults) <- list(rownames(results), cnz)
        newresults[,cnr] <- results
        results <- newresults
        for(i in 2:ncol(results))
          if(is.factor(results[,i])) results[,i] <- as.numeric(results[,i])
      }
      else if(all(!is.na(wh.zr))) {
        newz <- data.frame(matrix(NA,nrow=nrow(z),ncol=ncol(results)))
        dimnames(newz) <- list(rownames(z), cnr)
        newz[,cnz] <- z
        z <- newz
        for(i in 2:ncol(z))
          if(is.factor(z[,i])) z[,i] <- as.numeric(z[,i])
      }
      else {
        newnames <- c(cnr, cnz[is.na(wh.zr)])

        newresults <- data.frame(matrix(NA,nrow=nrow(results),ncol=length(newnames)))
        dimnames(newresults) <- list(rownames(results), newnames)
        newresults[,cnr] <- results
        results <- newresults
        for(i in 2:ncol(results))
          if(is.factor(results[,i])) results[,i] <- as.numeric(results[,i])
        
        newz <- data.frame(matrix(NA,nrow=nrow(z),ncol=length(newnames)))
        dimnames(newz) <- list(rownames(z), newnames)
        newz[,cnz] <- z
        z <- newz
        for(i in 2:ncol(z))
          if(is.factor(z[,i])) z[,i] <- as.numeric(z[,i])
      }
    }

    results <- rbind(results, z)
  }
  
  class(results) <- c("scanone","data.frame")
  attr(results,"method") <- method
  attr(results,"type") <- class(cross)[1]
  attr(results,"model") <- "twopart"
  results
}
