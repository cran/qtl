######################################################################
#
# scanone.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# last modified Dec, 2001
# first written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Hao Wu (The Jackson Lab) wrote the imputation method
#
# Part of the R/qtl package
# Contains: scanone, plot.scanone, scanone.perm
#           summary.scanone, print.summary.scanone
#
######################################################################

######################################################################
#
# scanone: scan genome, calculating LOD scores with single QTL model
#          (covariates are not allowed for models other than "normal")
#
######################################################################

scanone <-
function(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
         method=c("em","imp","hk","mr"), addcov=NULL, intcov=NULL,
         upper=FALSE, ties.random=FALSE, start=NULL, maxit=4000, tol=1e-4,
         n.perm, trace=TRUE)
{
  if(method=="im") # warning in case old terminology is used
    warning("Method \"im\" is now called \"em\"; running method \"imp\".")
  model <- match.arg(model)
  method <- match.arg(method)

  if(!missing(chr)) cross <- subset(cross, chr)

  # check phenotypes and covariates; drop individuals with missing values
  # in case of permutation test, only do checks once
  if(missing(n.perm) || n.perm > 0) {
    temp <- checkcovar(cross, pheno.col, addcov, intcov)
    cross <- temp[[1]]
    pheno <- temp[[2]]
    addcov <- temp[[3]]
    intcov <- temp[[4]]
    n.addcov <- temp[[5]]
    n.intcov <- temp[[6]]
  }
  else {
    pheno <- cross$pheno[,pheno.col]
    if(is.null(addcov)) n.addcov <- 0
    else n.addcov <- ncol(addcov)
    if(is.null(intcov)) n.intcov <- 0
    else n.intcov <- ncol(intcov)
  }
  n.chr <- nchr(cross)
  n.ind <- nind(cross)
  type <- class(cross)[1]

  # if n.perm specified, do a permutation test
  if(!missing(n.perm) && n.perm>0) {
    return(scanone.perm(cross, pheno.col, model, method, addcov,
                        intcov, upper, ties.random, start, maxit, 
                        tol, n.perm, trace))
  }

  if(model=="binary") {
    if(n.addcov > 0 || n.intcov > 0)
      warning("Covariates ignored for the binary model.")
    if(method=="imp" || method=="hk") {
      warning("Methods imp and hk not available for binary model; using em")
      method <- "em"
    }
    return(discan(cross,pheno.col, method, maxit, tol))
  }
  else if(model=="2part") {
    if(n.addcov > 0 || n.intcov > 0)
      warning("Covariates ignored for the two-part model.")
    if(method!="em") {
      warning("Only em method is available for the two-part model")
      method <- "em"
    }
    return(vbscan(cross, pheno.col, upper, method, maxit, tol))
  }
  else if(model=="np") {
    if(n.addcov > 0 || n.intcov > 0)
      warning("Covariates ignored for non-parametric interval mapping..")
    if(method!="em") {
      warning("Method argument ignored for non-parametric interval mapping.")
      method <- "em"
    }
  }
    
  # if non-parametric, convert phenotypes to ranks
  if(model=="np") {
    if(ties.random) {
      y <- pheno[!is.na(pheno)]
      y <- rank(y+runif(length(y))/(sd(y)*10^8))
      pheno[!is.na(pheno)] <- y
      correct <- 1
    }
    else {
      ties <- table(pheno)
      if(any(ties > 1)) {
        ties <- ties[ties>1]
        correct <- 1-sum(ties^3-ties)/(n.ind^3-n.ind)
      }
      else correct <- 1
      pheno <- rank(pheno)
    }
  }

  results <- NULL

  # starting points for interval mapping
  if(method=="em" && model=="normal") {
    if(is.null(start)) std.start <- 1
    else if(length(start)==1) std.start <- -1
    else std.start <- 0
  }

  # scan genome one chromosome at a time
  for(i in 1:n.chr) {

    # which type of cross is this?
    if(type == "f2") {
      if(class(cross$geno[[i]]) == "A") { # autosomal
        n.gen <- 3
        gen.names <- c("A","H","B")
      }
      else {                             # X chromsome 
        n.gen <- 2
        gen.names <- c("A","H","B") 
      }
    }
    else if(type == "bc") {
      n.gen <- 2
      gen.names <- c("A","H")
    }
    else if(type == "4way") {
      n.gen <- 4
      gen.names <- c("AC","AD","BC","BD")
    }
    else stop(paste("scanone not available for cross", type))

    # starting values for interval mapping
    if(method=="em" && model=="normal") {
      this.start <- rep(0,n.gen+1)
      if(std.start == 0) {
        if(length(start) < n.gen+1)
          stop(paste("Length of start argument should be 0, 1 or", n.gen+1))
        this.start <- c(start[1:n.gen],start[length(start)])
      }
    }

    # pull out reconstructed genotypes (mr)
    # or genotype probabilities (em or hk)
    if(method == "mr") {
      newgeno <- cross$geno[[i]]$data
      newgeno[is.na(newgeno)] <- 0 

      # discard partially informative genotypes
      if(type=="f2" || type=="f2ss") newgeno[newgeno>3] <- 0
      if(type=="4way") newgeno[newgeno>4] <- 0

      n.pos <- ncol(newgeno)
      map <- cross$geno[[i]]$map
      if(is.matrix(map)) {
        marnam <- colnames(map)
        map <- map[1,]
      }
      else marnam <- names(map)
    }
    else if(method == "imp") {
      if(is.na(match("draws",names(cross$geno[[i]])))) {
        # need to run sim.geno
        warning("First running sim.geno.")
        cross <- sim.geno(cross)
      }

      draws <- cross$geno[[i]]$draws
      n.pos <- ncol(draws)
      n.draws <- dim(draws)[3]

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$draws,"step"),
                        attr(cross$geno[[i]]$draws,"off.end"))
      if(is.matrix(map)) {
        marnam <- colnames(map)
        map <- map[1,]
      }
      else marnam <- names(map)
    }
    else {
      if(is.na(match("prob",names(cross$geno[[i]])))) {
        # need to run calc.genoprob
        warning("First running calc.genoprob.")
        cross <- calc.genoprob(cross)
      }
      genoprob <- cross$geno[[i]]$prob
      n.pos <- ncol(genoprob)

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$prob,"step"),
                        attr(cross$geno[[i]]$prob,"off.end"))
      if(is.matrix(map)) {
        marnam <- colnames(map)
        map <- map[1,]
      }
      else marnam <- names(map)
    }

    if(i==1 & (method=="hk" || method=="em")) { # get null log10 likelihood
      if(n.addcov > 0) resid0 <- lm(pheno ~ addcov)$resid
      else resid0 <- pheno - mean(pheno)
      if(method=="hk") nllik0 <- (n.ind/2)*log10(sum(resid0^2))
      else {
        sig0 <- sqrt(sum(resid0^2)/n.ind)
        nllik0 <- -sum(dnorm(resid0,0,sig0,log=TRUE))/log(10)
      }
    }

    # call the C function
    if(method == "mr") {
      z <- .C("R_scanone_mr",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.integer(newgeno),       # genotype data
              as.double(addcov),         # additive covariates
              as.integer(n.addcov),
              as.double(intcov),         # interactive covariates
              as.integer(n.intcov),
              as.double(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+2))),
              PACKAGE="qtl")
    }
    else if(method=="imp") {
      z <- .C("R_scanone_imp",
              as.integer(n.ind),
              as.integer(n.pos),
              as.integer(n.gen),
              as.integer(n.draws),
              as.integer(draws),
              as.double(addcov),
              as.integer(n.addcov),
              as.double(intcov),
              as.integer(n.intcov),
              as.double(pheno),
              result=as.double(rep(0,n.pos)),
              as.integer(1), # trim (for debugging purposes)
              as.integer(0), # direct (for debugging purposes)
              PACKAGE="qtl")
    }
    else if(method=="hk") { # Haley-Knott regression
      z <- .C("R_scanone_hk",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(addcov),         # additive covariates
              as.integer(n.addcov),
              as.double(intcov),         # interactive covariates
              as.integer(n.intcov), 
              as.double(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+2))),
              PACKAGE="qtl")
    }
    else if(method=="em" && model=="normal") { # interval mapping
      z <- .C("R_scanone_em",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(addcov),
              as.integer(n.addcov),
              as.double(intcov),
              as.integer(n.intcov),
              as.double(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+2))),
              as.integer(std.start),
              as.double(this.start),
              as.integer(maxit),
              as.double(tol),
              as.integer(0), # debugging trace off 
              PACKAGE="qtl")
    }
    else if(model=="np") { # non-parametric interval mapping
      z <- .C("R_scanone_np",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(pheno) ,         # phenotype data
              result=as.double(rep(0,n.pos)),
              PACKAGE="qtl")
    }
    else 
      stop(paste("Model", model, "with method", method, "not available"))

    z <- matrix(z$result,nrow=n.pos)

    # interval mapping without covariates:
    #   rescale log likelihood
    if(method!="imp" && n.addcov > 0)
      z <- z[,1,drop=FALSE]
    if(model == "np" && !ties.random)
      z <- z/correct  # correct for ties

    if(n.addcov==0 && n.intcov==0 && method != "imp"
       && model != "np") {
      # for the above cases, est'd coefficients are not returned
      if(type=="f2" && class(cross$geno[[i]])=="X") # add BB column
        z <- cbind(z[,1:3],rep(NA,n.pos),z[,4])
      colnames(z) <- c("lod",gen.names,"sigma")
    }
    else colnames(z) <- c("lod")
      
    w <- marnam
    o <- grep("^loc\-*[0-9]+",w)
    if(length(o) > 0) 
      w[o] <- paste(w[o],names(cross$geno)[i],sep=".c")
    rownames(z) <- w
    
    z <- as.data.frame(z)
    z <- cbind(chr=rep(names(cross$geno)[i],length(map)),
               pos=as.numeric(map), z)
    rownames(z) <- w
    results <- rbind(results,z)
  }

  # re-scale with null log10 likel for methods em and hk
  if((method=="em" && model=="normal") || method=="hk") 
    results[,3] <- nllik0 - results[,3]

  # replace any lod = NaN with 0
  results[is.na(results[,3]),3] <- 0

  class(results) <- c("scanone",class(results))
  attr(results,"method") <- method
  attr(results,"type") <- type
  attr(results,"model") <- model
  results
}


######################################################################
#
# plot.scanone: plot output from scanone
#
######################################################################

plot.scanone <- 
function(x,x2,x3,chr,incl.markers=TRUE,ylim,
         lty=c(1,2,3),col="black",lwd=2,add=FALSE,gap=25,
         main,...)
{
  if(length(dim(x))!=2)
    stop("Argument x must be a matrix or data.frame.")
  if(!missing(x2) && length(dim(x2))!=2)
    stop("Argument x2 must be a matrix or data.frame.")
  if(!missing(x3) && length(dim(x3))!=2)
    stop("Argument x3 must be matrices or data.frame.")

  second <- third <- TRUE
  if(missing(x2) && missing(x3)) 
     second <- third <- FALSE
  if(missing(x3))
    third <- FALSE
  if(missing(x2))
    second <- FALSE

  # rename things
  out <- x
  if(second) out2 <- x2
  if(third) out3 <- x3

  if(length(lty)==1) lty <- rep(lty,3)
  if(length(lwd)==1) lwd <- rep(lwd,3)
  if(length(col)==1) col <- rep(col,3)

  # pull out desired chromosomes
  if(missing(chr))
    chr <- unique(as.character(out[,1]))

  if(length(chr) == 0) chr <- sort(unique(out[,1]))
  else if(all(chr < 0)) { 
    a <- sort(unique(out[,1]))
    chr <- a[-match(-chr,a)]
  }
  out <- out[!is.na(match(out[,1],chr)),]
  if(second) out2 <- out2[!is.na(match(out2[,1],chr)),]
  if(third) out3 <- out3[!is.na(match(out3[,1],chr)),]
  
  # beginning and end of chromosomes
  temp <- grep("^loc\-*[0-9]+",rownames(out))
  if(length(temp)==0) temp <- out
  else temp <- out[-temp,]
  begend <- matrix(unlist(tapply(temp[,2],temp[,1],range)),ncol=2,byrow=TRUE)
  len <- begend[,2]-begend[,1]

  # locations to plot start of each chromosome
  start <- gap/2+c(0,cumsum(len+gap))-c(begend[,1],0)

  maxx <- sum(len+gap)
  maxy <- max(out[,3])
  if(second) maxy <- max(c(maxy,out2[,3]))
  if(third) maxy <- max(c(maxy,out3[,3]))

  # graphics parameters
  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=TRUE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))

  # make frame of plot
  if(missing(ylim)) ylim <- c(0,maxy)

  if(!add) {
    plot(0,0,ylim=ylim,xlim=c(0,maxx),type="n",
         xlab="Map position (cM)",ylab=dimnames(out)[[2]][3],
         ...)
    if(!missing(main)) title(main=main)
  }

  for(i in 1:length(chr)) {
    # plot first out
    x <- out[out[,1]==chr[i],2]+start[i]
    y <- out[out[,1]==chr[i],3]
    if(length(x)==1) {
      g <- max(gap/10,2)
      x <- c(x-g,x,x+g)
      y <- rep(y,3)
    }
    lines(x,y,lwd=lwd[1],lty=lty[1],col=col[1])

    # plot chromosome number
    a <- par("usr")
    if(!add) {
      tloc <- mean(c(min(x),max(x)))
      text(tloc,a[4]+(a[4]-a[3])*0.03,as.character(chr[i]))
      lines(rep(tloc,2),c(a[4],a[4]+(a[4]-a[3])*0.015))
    }

    # plot second out
    if(second) {
      x <- out2[out2[,1]==chr[i],2]+start[i]
      y <- out2[out2[,1]==chr[i],3]
      if(length(x)==1) {
        g <- max(gap/10,2)
        x <- c(x-g,x,x+g)
        y <- rep(y,3)
      }
      lines(x,y,lty=lty[2],col=col[2],lwd=lwd[2])
    }

    if(third) {
      x <- out3[out3[,1]==chr[i],2]+start[i]
      y <- out3[out3[,1]==chr[i],3]
      if(length(x)==1) {
        g <- max(gap/10,2)
        x <- c(x-g,x,x+g)
        y <- rep(y,3)
      }
      lines(x,y,lty=lty[3],col=col[3],lwd=lwd[3])
    }

    # plot lines at marker positions
    if(incl.markers && !add) {
      nam <- dimnames(out)[[1]][out[,1]==chr[i]]
      wh.genoprob <- (seq(along=nam))[grep("^loc\-*[0-9]+",nam)]
      if(length(wh.genoprob)==0) wh.genoprob <- seq(along=nam)
      else wh.genoprob <- (seq(along=nam))[-wh.genoprob]
      pos <- out[out[,1]==chr[i],2][wh.genoprob]+start[i]
      for(j in pos) 
	lines(c(j,j),c(a[3],a[3]+(a[4]-a[3])*0.02))
    }

  }

}

######################################################################
#
# scanone.perm: Permutation test of scanone
#
######################################################################

scanone.perm <-
function(cross, pheno.col=1, model=c("normal","binary","2part","np"),
         method=c("em","imp","hk","mr"), addcov=NULL, intcov=NULL,
         upper=FALSE, ties.random=FALSE, start=NULL, maxit=4000,
         tol=1e-4, n.perm=1000, trace=TRUE)
{
  if(method=="im") # warning in case old terminology is used
    warning("Method \"im\" is now called \"em\"; running method \"imp\".")
  method <- match.arg(method)
  model <- match.arg(model)

  if(model!="normal" && (!is.null(addcov) || !is.null(intcov))) {
    warning("Use of covariates not available for method np")
    addcov <- intcov <- NULL
  }

  n.ind <- nind(cross)

  addcovp <- intcovp <- NULL
  if(!is.null(addcov)) addcov <- as.matrix(addcov)
  if(!is.null(intcov)) intcov <- as.matrix(intcov)

  if(model=="2part") res <- matrix(ncol=3,nrow=n.perm)
  else res <- 1:n.perm

  if(trace) { # if trace, print out a tracing information
    # rnd: how often to print tracing information
    if(trace > 1) rnd <- 1
    else {
      if(n.perm >= 1000) rnd <- 20
      else if(n.perm >= 100) rnd <- 5
      else rnd <- 1
    }
  }

  for(i in 1:n.perm) {
    if(trace && i/rnd == round(i/rnd))
      cat("Permutation", i, "\n")

    o <- sample(1:n.ind)
    cross$pheno <- cross$pheno[o,,drop=FALSE]
    if(!is.null(addcov)) addcovp <- addcov[o,,drop=FALSE]
    if(!is.null(intcov)) intcovp <- intcov[o,,drop=FALSE]
    tem <- scanone(cross,,pheno.col,model,method,addcovp,
                   intcovp,upper,ties.random,start,maxit,tol,
                   n.perm= -1)
    if(model=="2part") res[i,] <- apply(tem[,3:5],2,max,na.rm=TRUE)
    else res[i] <- max(tem[,3],na.rm=TRUE)
  }

  if(model=="2part") colnames(res) <- c("LOD.p.mu", "LOD.p", "LOD.mu")

  attr(res,"method") <- method
  attr(res,"model") <- model
  attr(res,"type") <- class(cross)[1]
  res
}

# give, for each chromosome, the position with the maximum LOD
summary.scanone <-
function(object,threshold=0,...)
{
  out <- lapply(split(object,object[,1]),
                   function(b) b[b[,3]==max(b[,3]),])
  results <- out[[1]]
  if(length(out) > 1)
    for(i in 2:length(out))
      results <- rbind(results,out[[i]])
  class(results) <- c("summary.scanone","data.frame")
  if(!any(results[,3] >= threshold)) {
    cat("    There were no LOD peaks above the threshold.\n")
    invisible()
  }
  else return(results[results[,3] >= threshold,])
}

# print output of summary.scanone
print.summary.scanone <-
function(x,...)
{
  x[,-(1:2)] <- round(data.frame(x[,-(1:2)]),6)
  cat("\n")
  print.data.frame(x,digits=2)
  cat("\n")
}

# end of scanone.R
