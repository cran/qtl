######################################################################
#
# discan.R
#
# copyright (c) 2001-4, Karl W Broman, Johns Hopkins University
# last modified Aug, 2004
# first written Oct, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: discan
#
######################################################################

######################################################################
#
# discan: scan genome, calculating LOD scores with single QTL model
#         for a dichotomous trait 
#
######################################################################

discan <-
function(cross, pheno.col=1, method=c("em","mr"),
         maxit=4000, tol=1e-4)
{
  method <- match.arg(method)

  # check phenotypes
  if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
  if(pheno.col < 1 || pheno.col > nphe(cross))
    stop("Specified phenotype column is invalid.")

  pheno <- cross$pheno[,pheno.col]
  keep.ind <- (1:length(pheno))[!is.na(pheno)]
  pheno <- pheno[keep.ind]
  n.ind <- length(keep.ind)
  n.chr <- nchr(cross)
  type <- class(cross)[1]

  u <- unique(pheno)
  if(any(u != 0 && u != 1))
    stop("Phenotypes must be either 0 or 1.")

  if(method == "em") {
    p <- mean(pheno)
    n1 <- sum(pheno==1)
    n0 <- sum(pheno==0)
    if(n1==0 || n0==0) llik0 <- 0
    else llik0 <- n1*log10(p) + n0*log10(1-p)
  }

  results <- NULL

  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {

    chrtype <- class(cross$geno[[i]])
    if(chrtype=="X") sexpgm <- getsex(cross)
    else sexpgm <- NULL

    # get genotype names
    gen.names <- getgenonames(type,chrtype,"full",sexpgm)
    n.gen <- length(gen.names)

    # pull out genotype data (mr)
    # or genotype probabilities (em)
    if(method == "mr") {
      cfunc <- "R_discan_mr"
      newgeno <- cross$geno[[i]]$data
      newgeno <- newgeno[keep.ind,]
      newgeno[is.na(newgeno)] <- 0 

      # discard partially informative genotypes
      if(type=="f2" || type=="f2ss") newgeno[newgeno>3] <- 0
      if(type=="4way") newgeno[newgeno>4] <- 0

      # revise X chromosome genotypes
      if(chrtype=="X" && (type=="bc" || type=="f2" || type=="f2ss"))
         newgeno <- reviseXdata(type, "full", sexpgm, geno=newgeno)

      n.pos <- ncol(newgeno)
      map <- cross$geno[[i]]$map
      if(is.matrix(map)) map <- map[1,]
    }
    else {
      if(is.na(match("prob",names(cross$geno[[i]])))) { # need to run calc.genoprob
        warning("First running calc.genoprob.")
        cross <- calc.genoprob(cross)
      }
      genoprob <- cross$geno[[i]]$prob
      n.pos <- ncol(genoprob)
      genoprob <- genoprob[keep.ind,,]

      # revise X chromosome genotypes
      if(chrtype=="X" && (type=="bc" || type=="f2" || type=="f2ss"))
         genoprob <- reviseXdata(type, "full", sexpgm, prob=genoprob)

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$prob,"step"),
                        attr(cross$geno[[i]]$prob,"off.end"))
      if(is.matrix(map)) map <- map[1,]

      cfunc <- "R_discan_im"
    }

    # call the C function
    if(method == "mr") 
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.integer(newgeno),       # genotype data
              as.integer(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+1))),
              PACKAGE="qtl")

    else  # interval mapping
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.integer(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+1))),
              as.integer(maxit),
              as.double(tol),
              PACKAGE="qtl")
    z <- matrix(z$result,nrow=n.pos)

    if(method == "em") z[,1] <- z[,1] - llik0
    z[is.na(z[,1]),1] <- 0
    colnames(z) <- c("lod",gen.names)
      
    w <- names(map)
    o <- grep("^loc\-*[0-9]+",w)
    if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
      w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
    rownames(z) <- w
    
    z <- as.data.frame(z)
    z <- cbind(chr=rep(names(cross$geno)[i],length(map)), pos=map, z)
    rownames(z) <- w


    # get null log10 likelihood for the X chromosome
    if(chrtype=="X") {

      # determine which covariates belong in null hypothesis
      temp <- scanoneXnull(type, sexpgm)
      adjustX <- temp$adjustX
      dfX <- temp$dfX
      sexpgmcovar <- temp$sexpgmcovar
      sexpgmcovar.alt <- temp$sexpgmcovar.alt      

      if(adjustX) { # get LOD-score adjustment
        n.gen <- ncol(sexpgmcovar)+1

        nullz <- .C("R_discan_mr",
            as.integer(n.ind),
            as.integer(1),
            as.integer(n.gen),
            as.integer(sexpgmcovar.alt),
            as.integer(pheno),
            result=as.double(rep(0,n.gen+1)),
            PACKAGE="qtl")

        # adjust LOD curve
        z[,3] <- z[,3] - nullz$result[1]
      }
    } 

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

  # sort the later columns
  neworder <- c(colnames(results)[1:3],sort(colnames(results)[-(1:3)]))
  results <- results[,neworder]

  class(results) <- c("scanone","data.frame")
  attr(results,"method") <- method
  attr(results,"type") <- type
  attr(results,"model") <- "binary"
  results
}

# end of discan.R
