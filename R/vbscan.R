######################################################################
#
# vbscan.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# July, 2001; May, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: vbscan, vbscan.perm
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
function(cross, chr, pheno.col=1, upper=FALSE,
	 maxit=1000, tol=1e-8)
{
  if(!missing(chr)) cross <- pull.chr(cross, chr)

  # check arguments are okay
  if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
  if(pheno.col > nphe(cross))
    stop("Specified phenotype column exceeds the number of phenotypes")
  y <- cross$pheno[,pheno.col]

  # remove individuals with missing phenotypes
  if(any(is.na(y))) {
    drop <- (1:length(y))[is.na(y)]
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

  results <- NULL
  for(i in 1:length(cross$geno)) {
    # make sure inferred genotypes or genotype probabilities are available
    if(is.na(match("prob",names(cross$geno[[i]])))) {
      cat(" -Calculating genotype probabilities\n")
      cross <- calc.genoprob(cross)
    }

    n.pos <- dim(cross$geno[[i]]$prob)[2]
    n.ind <- length(y)
    n.gen <- dim(cross$geno[[i]]$prob)[3]

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
                       paste("pi",c("AA","AB","BB")[1:n.gen]),
                       paste("mu",c("AA","AB","BB")[1:n.gen]),
                       "sigma")

    # the following is for the case where there is a sex chromosome
    if(!is.null(results) && ncol(results) != ncol(res)) {
      if(ncol(results) > ncol(res)) {
        o <- match(colnames(res),colnames(results))
        missing <- (1:ncol(results))[-o]
        temp <- as.data.frame(matrix(ncol=ncol(results),nrow=nrow(res)))
        colnames(temp) <- colnames(results)
        rownames(temp) <- rownames(res)
        temp[,o] <- res
        res <- temp
      }
      else {
        o <- match(colnames(results),colnames(res))
        missing <- (1:ncol(res))[-o]
        temp <- as.data.frame(matrix(ncol=ncol(res),nrow=nrow(results)))
        colnames(temp) <- colnames(res)
        rownames(temp) <- rownames(results)
        temp[,o] <- results
        results <- temp
      }
    }
    
    results <- rbind(results,res)
  }
  
  class(results) <- c("scanone","data.frame")
  results
}


# permutation test for vbscan
vbscan.perm <-
function(cross, chr, pheno.col=1, upper=FALSE,
	 n.perm=1000, maxit=1000, tol=1e-8)
{
  if(!missing(chr)) cross <- pull.chr(cross, chr)

  n.ind <- nind(cross)

  res <- matrix(ncol=3,nrow=n.perm)
  for(i in 1:n.perm) {
    cross$pheno <- as.matrix(cross$pheno[sample(1:n.ind),])
    tem <- vbscan(cross,pheno.col=pheno.col,upper=upper,
                      maxit=maxit,tol=tol)
#    print(c(i,dim(tem)))
    res[i,] <- apply(tem[,3:5],2,max)
  }
  colnames(res) <- c("lod","lod.p","lod.mu")
  res
}
