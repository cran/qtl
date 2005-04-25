######################################################################
#
# discan.R
#
# copyright (c) 2001-5, Karl W Broman, Johns Hopkins University
# last modified Apr, 2005
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
function(cross, pheno, method=c("em","mr"),
         addcovar=NULL, intcovar=NULL, maxit=4000, tol=1e-4,
         verbose=FALSE)
{
  method <- match.arg(method)

  n.ind <- nind(cross)
  n.chr <- nchr(cross)
  type <- class(cross)[1]
  if(is.null(addcovar)) {
    n.addcovar <- 0
    addcovar <- 0
  }
  else n.addcovar <- ncol(addcovar)
  if(is.null(intcovar)) {
    n.intcovar <- 0
    intcovar <- 0
  }
  else n.intcovar <- ncol(intcovar)

  if(method=="mr" && n.addcovar+n.intcovar>0)  {
    warning("Covariates ignored with method=\"mr\"; use \"em\" instead")
    n.addcovar <- n.intcovar <- addcovar <- intcovar <- 0
  }

  u <- unique(pheno)
  if(any(u != 0 & u != 1))
    stop("Phenotypes must be either 0 or 1.")

  # get null log liklihood
  if(n.addcovar > 0)
    nullfit <- glm(pheno ~ addcovar, family=binomial(link=logit))
  else
    nullfit <- glm(pheno ~ 1, family=binomial(link=logit))
  fitted <- nullfit$fitted
  nullcoef <- nullfit$coef
  llik0 <- sum(pheno*log10(fitted) + (1-pheno)*log10(1-fitted))

  results <- NULL

  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {

    chrtype <- class(cross$geno[[i]])
    if(chrtype=="X") {
      sexpgm <- getsex(cross)
      ac <- revisecovar(sexpgm,addcovar)
      n.ac <- ifelse(is.null(ac),0,ncol(ac))
      ic <- revisecovar(sexpgm,intcovar)
      n.ic <- ifelse(is.null(ic),0,ncol(ic))
    }
    else {
      sexpgm <- NULL
      ac <- addcovar
      n.ac <- n.addcovar
      ic <- intcovar
      n.ic <- n.intcovar
    }

    # get genotype names
    gen.names <- getgenonames(type,chrtype,"full",sexpgm)
    n.gen <- length(gen.names)

    # pull out genotype data (mr)
    # or genotype probabilities (em)
    if(method == "mr") {
      newgeno <- cross$geno[[i]]$data
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

      z <- .C("R_discan_mr",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.integer(newgeno),       # genotype data
              as.integer(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+1))),
              PACKAGE="qtl")

    }
    else {
      if(is.na(match("prob",names(cross$geno[[i]])))) { # need to run calc.genoprob
        warning("First running calc.genoprob.")
        cross <- calc.genoprob(cross)
      }
      genoprob <- cross$geno[[i]]$prob
      n.pos <- ncol(genoprob)

      # revise X chromosome genotypes
      if(chrtype=="X" && (type=="bc" || type=="f2" || type=="f2ss"))
         genoprob <- reviseXdata(type, "full", sexpgm, prob=genoprob)

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$prob,"step"),
                        attr(cross$geno[[i]]$prob,"off.end"))
      if(is.matrix(map)) map <- map[1,]

      if(n.ac + n.ic > 0) {

        start <- rep(nullcoef[1],n.gen)
        if(n.ac > 0)
          start <- c(start, nullcoef[-1])
        if(n.ic > 0)
          start <- c(start, rep(0, n.ic*(n.gen-1)))

        z <- .C("R_discan_covar",
                as.integer(n.ind),         # number of individuals
                as.integer(n.pos),         # number of markers
                as.integer(n.gen),         # number of possible genotypes
                as.double(genoprob),       # genotype probabilities
                as.double(ac),
                as.integer(n.ac),
                as.double(ic),
                as.integer(n.ic),
                as.integer(pheno),          # phenotype data
                as.double(start),
                result=as.double(rep(0,n.pos)),
                as.integer(maxit),
                as.double(tol),
                as.integer(verbose),
                PACKAGE="qtl")
      }
      else {
        z <- .C("R_discan_im",
                as.integer(n.ind),         # number of individuals
                as.integer(n.pos),         # number of markers
                as.integer(n.gen),         # number of possible genotypes
                as.double(genoprob),       # genotype probabilities
                as.integer(pheno),          # phenotype data
                result=as.double(rep(0,n.pos*(n.gen+1))),
                as.integer(maxit),
                as.double(tol),
                PACKAGE="qtl")
      }

    }
    z <- matrix(z$result,nrow=n.pos)

    if(method == "em") z[,1] <- z[,1] - llik0
    z[is.na(z[,1]),1] <- 0

    if(n.ac + n.ic > 0)
      colnames(z) <- c("lod")
    else
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
    adjustX <- FALSE
    if(chrtype=="X") {

      # determine which covariates belong in null hypothesis
      temp <- scanoneXnull(type, sexpgm)
      adjustX <- temp$adjustX
      dfX <- temp$dfX
      sexpgmcovar <- temp$sexpgmcovar
      sexpgmcovar.alt <- temp$sexpgmcovar.alt      

      if(adjustX) { # get LOD-score adjustment
        if(n.ac > 0) 
          nullfit <- glm(pheno ~ ac+sexpgmcovar,
                         family=binomial(link=logit))
        else 
          nullfit <- glm(pheno ~ sexpgmcovar,
                         family=binomial(link=logit))
        fitted <- nullfit$fitted
        llik0X <- sum(pheno*log10(fitted) + (1-pheno)*log10(1-fitted))

        # adjust LOD curve
        z[,3] <- z[,3] - (llik0X - llik0)
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
  } # loop over chromosomes

  if(ncol(results) > 3) {
    # sort the later columns
    neworder <- c(colnames(results)[1:3],sort(colnames(results)[-(1:3)]))
    results <- results[,neworder]
  }

  class(results) <- c("scanone","data.frame")
  attr(results,"method") <- method
  attr(results,"type") <- type
  attr(results,"model") <- "binary"
  attr(results,"null.log10.lik") <- llik0
  if(adjustX) 
    attr(results,"null.log10.lik.X") <- llik0X
  results
}

# end of discan.R
