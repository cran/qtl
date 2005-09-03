#####################################################################
#
# scanone.R
#
# copyright (c) 2001-5, Karl W Broman, Johns Hopkins University
# last modified Sep, 2005
# first written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Hao Wu (The Jackson Lab) wrote the imputation method
#
# Part of the R/qtl package
# Contains: scanone, plot.scanone, scanone.perm,
#           summary.scanone, print.summary.scanone,
#           max.scanone
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
         method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
         addcovar=NULL, intcovar=NULL, weights=NULL,
         upper=FALSE, ties.random=FALSE,
         start=NULL, maxit=4000, tol=1e-4, n.perm, verbose)
{
  model <- match.arg(model)
  method <- match.arg(method)

  if(missing(verbose)) {
    if(!missing(n.perm) && n.perm > 0) verbose <- TRUE
    else verbose <- FALSE
  }
    
  if(!missing(chr)) cross <- subset(cross, chr)
  if(missing(n.perm)) n.perm <- 0

  # check phenotypes and covariates; drop individuals with missing values
  # in case of permutation test, only do checks once
  if(n.perm >= 0) {
    temp <- checkcovar(cross, pheno.col, addcovar, intcovar)
    cross <- temp[[1]]
    pheno <- temp[[2]]
    addcovar <- temp[[3]]
    intcovar <- temp[[4]]
    n.addcovar <- temp[[5]]
    n.intcovar <- temp[[6]]
  }
  else {
    pheno <- cross$pheno[,pheno.col]
    if(is.null(addcovar)) n.addcovar <- 0
    else n.addcovar <- ncol(addcovar)
    if(is.null(intcovar)) n.intcovar <- 0
    else n.intcovar <- ncol(intcovar)
  }
  n.chr <- nchr(cross)
  n.ind <- nind(cross)
  type <- class(cross)[1]

  # if n.perm specified, do a permutation test
  if(n.perm>0) {
    return(scanone.perm(cross, pheno.col, model, method, addcovar,
                        intcovar, weights, upper, ties.random,
                        start, maxit, tol, n.perm, verbose))
  }

  # fill in missing genotypes with imputed values
  if(n.perm==0) { # not in the midst of permutations
    if(method=="mr-argmax")
      cross <- fill.geno(cross,method="argmax")
    if(method=="mr-imp")
      cross <- fill.geno(cross,method="imp")
  }

  # weights for model="normal"
  if(model != "normal") {
    if(!is.null(weights) && !all(weights==1)) 
      warning("weights used only for normal model.")
  }
  else {
    if(is.null(weights))
      weights <- rep(1, nind(cross))
    else 
      if(length(weights) != nind(cross))
        stop("weights should either be NULL or a vector of length n.ind")
    if(any(weights) <= 0)
      stop("weights should be entirely positive")
    weights <- sqrt(weights)
  }

  if(model=="binary") {
    if(method=="imp" || method=="hk") {
      warning("Methods imp and hk not available for binary model; using em")
      method <- "em"
    }
    return(discan(cross, pheno, method, addcovar, intcovar, maxit, tol, verbose))
  }
  else if(model=="2part") {
    if(n.addcovar > 0 || n.intcovar > 0)
      warning("Covariates ignored for the two-part model.")
    if(method!="em") {
      warning("Only em method is available for the two-part model")
      method <- "em"
    }
    return(vbscan(cross, pheno.col, upper, method, maxit, tol))
  }
  else if(model=="np") {
    if(n.addcovar > 0 || n.intcovar > 0)
      warning("Covariates ignored for non-parametric interval mapping.")
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

    # starting values for interval mapping
    if(method=="em" && model=="normal") {
      this.start <- rep(0,n.gen+1)
      if(std.start == 0) {
        if(length(start) < n.gen+1) 
          stop("Length of start argument should be 0, 1 or ", n.gen+1)
        this.start <- c(start[1:n.gen],start[length(start)])
      }
    }

    # pull out reconstructed genotypes (mr)
    # or imputations (imp)
    # or genotype probabilities (em or hk)
    if(method=="mr" || method=="mr-imp" || method=="mr-argmax") {
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

      # revise X chromosome genotypes
      if(chrtype=="X" && (type=="bc" || type=="f2" || type=="f2ss"))
         draws <- reviseXdata(type, "full", sexpgm, draws=draws)

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

      # revise X chromosome genotypes
      if(chrtype=="X" && (type=="bc" || type=="f2" || type=="f2ss"))
         genoprob <- reviseXdata(type, "full", sexpgm, prob=genoprob)

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$prob,"step"),
                        attr(cross$geno[[i]]$prob,"off.end"))
      if(is.matrix(map)) {
        marnam <- colnames(map)
        map <- map[1,]
      }
      else marnam <- names(map)
    }

    # call the C function
    if(method == "mr" || method=="mr-imp" || method=="mr-argmax") 
      z <- .C("R_scanone_mr",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.integer(newgeno),       # genotype data
              as.double(ac),       # additive covariates
              as.integer(n.ac),
              as.double(ic),       # interactive covariates
              as.integer(n.ic),
              as.double(pheno),          # phenotype data
              as.double(weights),        # weights
              result=as.double(rep(0,n.pos*(n.gen+2))),
              PACKAGE="qtl")

    else if(method=="imp") 
      z <- .C("R_scanone_imp",
              as.integer(n.ind),
              as.integer(n.pos),
              as.integer(n.gen),
              as.integer(n.draws),
              as.integer(draws),
              as.double(ac),
              as.integer(n.ac),
              as.double(ic),
              as.integer(n.ic),
              as.double(pheno),
              as.double(weights),
              result=as.double(rep(0,n.pos)),
              as.integer(1), # trim (for debugging purposes)
              as.integer(0), # direct (for debugging purposes)
              PACKAGE="qtl")
    
    else if(method=="hk")  # Haley-Knott regression
      z <- .C("R_scanone_hk",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(ac),         # additive covariates
              as.integer(n.ac),
              as.double(ic),         # interactive covariates
              as.integer(n.ic), 
              as.double(pheno),          # phenotype data
              as.double(weights),
              result=as.double(rep(0,n.pos*(n.gen+2))),
              PACKAGE="qtl")
   
    else if(method=="em" && model=="normal")  # interval mapping
      z <- .C("R_scanone_em",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(ac),
              as.integer(n.ac),
              as.double(ic),
              as.integer(n.ic),
              as.double(pheno),          # phenotype data
              as.double(weights),
              result=as.double(rep(0,n.pos*(n.gen+2))),
              as.integer(std.start),
              as.double(this.start),
              as.integer(maxit),
              as.double(tol),
              as.integer(0), # debugging verbose off 
              PACKAGE="qtl")

    else if(model=="np")  # non-parametric interval mapping
      z <- .C("R_scanone_np",
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(pheno) ,         # phenotype data
              result=as.double(rep(0,n.pos)),
              PACKAGE="qtl")
    
    else  {
      err <- paste("Model", model, "with method", method, "not available")
      stop(err)
    }

    z <- matrix(z$result,nrow=n.pos)

    # interval mapping without covariates:
    #   rescale log likelihood
    if(method!="imp" && n.ac > 0)
      z <- z[,1,drop=FALSE]
    if(model == "np" && !ties.random)
      z <- z/correct  # correct for ties

    if(n.ac==0 && n.ic==0 && method != "imp"
       && model != "np") 
      colnames(z) <- c("lod",gen.names,"sigma")
    else colnames(z) <- c("lod")
      
    w <- marnam
    o <- grep("^loc\-*[0-9]+",w)
    if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
      w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
    rownames(z) <- w
    
    z <- as.data.frame(z)
    z <- cbind(chr=rep(names(cross$geno)[i],length(map)),
               pos=as.numeric(map), z)
    rownames(z) <- w

    # get null log10 likelihood
    if(i==1 & model != "np") {
      if(n.ac > 0)
        resid0 <- lm(pheno ~ ac, weights=weights^2)$resid
      else 
        resid0 <- lm(pheno ~ 1, weights=weights^2)$resid
      if(method=="hk") nllik0 <- (n.ind/2)*log10(sum((resid0*weights)^2))
      else {
        sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
        nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
      }
    }

    # re-scale with null log10 likel for methods em and hk
    if((method=="em" && model=="normal") || method=="hk") 
      z[,3] <- nllik0 - z[,3]

    # get null log10 likelihood for the X chromosome
    if(chrtype=="X") {

      # determine which covariates belong in null hypothesis
      temp <- scanoneXnull(type, sexpgm)
      adjustX <- temp$adjustX
      dfX <- temp$dfX
      sexpgmcovar <- temp$sexpgmcovar
      
      if(adjustX) {
        if(model == "np") {
          sexpgmcovar <- factor(apply(sexpgmcovar,1,paste,collapse=":"))
          nllikX <- kruskal.test(pheno ~ sexpgmcovar)$stat/(2*log(10))
          z[,3] <- z[,3] - nllikX
        }
        else if(method=="mr") {
          for(s in 1:ncol(newgeno)) {
            wh <- newgeno[,s] != 0
            
            if(n.ac > 0) {
              residX <- lm(pheno ~ ac+sexpgmcovar, weights=weights^2,subset=wh)$resid
              resid0 <- lm(pheno ~ ac, weights=weights^2,subset=wh)$resid
            }
            else {
              residX <- lm(pheno ~ sexpgmcovar, weights=weights^2,subset=wh)$resid
              resid0 <- lm(pheno ~ 1, weights=weights^2,subset=wh)$resid
            }
            nllikX <- (sum(wh)/2)*log10(sum((residX*weights[wh])^2))
            nllik0 <- (sum(wh)/2)*log10(sum((resid0*weights[wh])^2))

            # rescale LOD score
            z[s,3] <- z[s,3] + nllikX - nllik0
          }
        }
        else {
          if(n.ac > 0) {
            outX <- lm(pheno ~ ac+sexpgmcovar, weights=weights^2)
            residX <- outX$resid
            # perhaps revise the dfX, if some columns got dropped
            dfX <- dfX - (ncol(sexpgmcovar)+n.ac - (outX$rank-1))
          }
          else 
            residX <- lm(pheno ~ sexpgmcovar, weights=weights^2)$resid

          if(method=="hk") nllikX <- (n.ind/2)*log10(sum((residX*weights)^2))
          else {
            if(method=="imp") {
              if(n.ac > 0) {
                out0 <- lm(pheno ~ ac, weights=weights^2)
                resid0 <- out0$resid
              }
              else {
                out0 <- lm(pheno ~ 1, weights=weights^2)
                resid0 <- out0$resid
              }
              
              sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
              nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
            }
            sigX <- sqrt(sum((residX*weights)^2)/n.ind)
            nllikX <- -sum(dnorm(residX,0,sigX/weights,log=TRUE))/log(10)
          }
          # rescale LOD score
          zz <- z
          z[,3] <- z[,3] + nllikX - nllik0
        }
      }
    }

    # replace missing or negative LODs with 0
    z[is.na(z[,3]) | z[,3]<0, 3] <- 0

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
        for(s in 2:ncol(results))
          if(is.factor(results[,s])) results[,s] <- as.numeric(results[,s])
      }
      else if(all(!is.na(wh.zr))) {
        newz <- data.frame(matrix(NA,nrow=nrow(z),ncol=ncol(results)))
        dimnames(newz) <- list(rownames(z), cnr)
        newz[,cnz] <- z
        z <- newz
        for(s in 2:ncol(z))
          if(is.factor(z[,s])) z[,s] <- as.numeric(z[,s])
      }
      else {
        newnames <- c(cnr, cnz[is.na(wh.zr)])

        newresults <- data.frame(matrix(NA,nrow=nrow(results),ncol=length(newnames)))
        dimnames(newresults) <- list(rownames(results), newnames)
        newresults[,cnr] <- results
        results <- newresults
        for(s in 2:ncol(results))
          if(is.factor(results[,s])) results[,s] <- as.numeric(results[,s])
        
        newz <- data.frame(matrix(NA,nrow=nrow(z),ncol=length(newnames)))
        dimnames(newz) <- list(rownames(z), newnames)
        newz[,cnz] <- z
        z <- newz
        for(s in 2:ncol(z))
          if(is.factor(z[,s])) z[,s] <- as.numeric(z[,s])
      }
    }

    results <- rbind(results,z)
  } # end loop over chromosomes

  if(n.addcovar + n.intcovar > 0)
    results <- results[,1:3]
  
  # sort the later columns
  if(ncol(results) > 3) {
    neworder <- c(colnames(results)[1:3],sort(colnames(results)[-(1:3)]))
    results <- results[,neworder]
  }

  class(results) <- c("scanone","data.frame")
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
function(x,x2,x3,chr,lodcolumn=3,incl.markers=TRUE,xlim, ylim,
         lty=1,col=c("black","blue","red"),lwd=2,add=FALSE,gap=25,
         main, mtick=c("line", "triangle"), ...)
{
  mtick <- match.arg(mtick)

  if(length(dim(x))!=2)
    stop("Argument x must be a matrix or data.frame.")
  if(!missing(x2) && length(dim(x2))!=2)
    stop("Argument x2 must be a matrix or data.frame.")
  if(!missing(x3) && length(dim(x3))!=2)
    stop("Argument x3 must be a matrix or data.frame.")

  if(length(lodcolumn)==1) 
    lodcolumn <- rep(lodcolumn,3)[1:3]
  else if(length(lodcolumn)==2) {
    if(missing(x2)) x2 <- x
    lodcolumn <- lodcolumn[c(1,2,3)]
  }
  else {
    if(missing(x2)) x2 <- x
    if(missing(x3)) x3 <- x
  }

  second <- third <- TRUE
  if(missing(x2) && missing(x3)) 
     second <- third <- FALSE
  if(missing(x3))
    third <- FALSE
  if(missing(x2))
    second <- FALSE

  # rename things and turn into data frames
  out <- x[,c(1:2,lodcolumn[1])]
  if(second) out2 <- x2[,c(1:2,lodcolumn[2])]
  if(third) out3 <- x3[,c(1:2,lodcolumn[3])]
  if(length(lty)==1) lty <- rep(lty,3)
  if(length(lwd)==1) lwd <- rep(lwd,3)
  if(length(col)==1) col <- rep(col,3)

  # pull out desired chromosomes
  if(missing(chr) || length(chr)==0) 
    chr <- unique(as.character(out[,1]))
  else if(all(chr < 0)) { 
    a <- sort(unique(out[,1]))
    chr <- a[-match(-chr,a)]
  }

  u <- is.na(match(chr,unique(out[,1])))
  if(all(u))
    stop("Chromosome(s) to plot were not matched to those in the scanone output.")
  else if(any(u)) {
    warning(paste("Chromosome(s)",chr[u],"were not found.",sep=" "))
    chr <- chr[!u]
  }

  out <- out[!is.na(match(out[,1],chr)),]
  if(second) out2 <- out2[!is.na(match(out2[,1],chr)),]
  if(third) out3 <- out3[!is.na(match(out3[,1],chr)),]
  
  onechr <- FALSE
  if(length(chr) == 1) {
    gap <- 0
    onechr <- TRUE 
  }

  # beginning and end of chromosomes
#  temp <- grep("^c[0-9A-Za-z]+\.loc\-*[0-9]+",rownames(out))
#  if(length(temp)==0) temp <- out
#  else temp <- out[-temp,]
  temp <- out
  begend <- matrix(unlist(tapply(temp[,2],temp[,1],range)),ncol=2,byrow=TRUE)
  rownames(begend) <- unique(out[,1])
  begend <- begend[as.character(chr),,drop=FALSE]
  len <- begend[,2]-begend[,1]

  # locations to plot start of each chromosome
  if(!onechr) 
    start <- c(0,cumsum(len+gap))-c(begend[,1],0)
  else start <- 0

  maxx <- sum(len+gap)-gap
  maxy <- max(out[,3],na.rm=TRUE)
  if(second) maxy <- max(c(maxy,out2[,3]),na.rm=TRUE)
  if(third) maxy <- max(c(maxy,out3[,3]),na.rm=TRUE)

  # graphics parameters
  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=FALSE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))

  # make frame of plot
  if(missing(ylim)) ylim <- c(0,maxy)
  if(missing(xlim)) {
    if(onechr) xlim <- c(0,max(out[,2]))
    else xlim <- c(0,maxx)
  }
  
  if(!add) {
    if(onechr) {
      plot(0,0,ylim=ylim,xlim=xlim,type="n",
           xlab="Map position (cM)",ylab=dimnames(out)[[2]][3],
           ...)
    }
    else {
      plot(0,0,ylim=ylim,xlim=xlim,type="n",xaxt="n",
           xlab="",ylab=dimnames(out)[[2]][3],
           ...)
    }
    if(!missing(main)) title(main=main)
  }

  # initialize xtick and xtickmark
  xtick <- NULL
  xticklabel <- NULL
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
#    a <- par("usr")
    if(!add && !onechr) {
      tloc <- mean(c(min(x),max(x)))
#      text(tloc,a[3]-(a[4]-a[3])*0.05,as.character(chr[i]))
#      lines(rep(tloc,2),c(a[3],a[3]-(a[4]-a[3])*0.015))
      xtick <- c(xtick, tloc)
      xticklabel <- c(xticklabel, as.character(chr[i]))
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

    # plot lines or triangles at marker positions
    if(incl.markers && !add) {
      nam <- dimnames(out)[[1]][out[,1]==chr[i]]
#      wh.genoprob <- (seq(along=nam))[grep("^loc\-*[0-9]+",nam)]
      wh.genoprob <- (seq(along=nam))[grep("^c[0-9A-Za-z]+\.loc\-*[0-9]+",nam)]
      if(length(wh.genoprob)==0) wh.genoprob <- seq(along=nam)
      else wh.genoprob <- (seq(along=nam))[-wh.genoprob]
      pos <- out[out[,1]==chr[i],2][wh.genoprob]+start[i]
      if(mtick=="line")
        rug(pos, 0.02, quiet=TRUE)
      else {
        a <- par("usr")
        points(pos, rep(a[3]+diff(a[3:4])*0.04, length(pos)), pch=17, cex=1.5)
      }
      #for(j in pos)
      #  lines(c(j,j),c(a[3],a[3]+(a[4]-a[3])*0.02))
    }
  }
  # draw the axis
  if(!add && !onechr) 
    axis(1, at=xtick, labels=xticklabel)

  invisible()

}

######################################################################
#
# scanone.perm: Permutation test of scanone
#
######################################################################

scanone.perm <-
function(cross, pheno.col=1, model=c("normal","binary","2part","np"),
         method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
         addcovar=NULL, intcovar=NULL, weights=NULL,
         upper=FALSE, ties.random=FALSE,
         start=NULL, maxit=4000, tol=1e-4, n.perm=1000, verbose=TRUE)
{
  method <- match.arg(method)
  model <- match.arg(model)

  if(model!="normal" && (!is.null(addcovar) || !is.null(intcovar))) {
    warning("Use of covariates not available for method np")
    addcovar <- intcovar <- NULL
  }

  n.ind <- nind(cross)

  addcovarp <- intcovarp <- NULL
  if(!is.null(addcovar)) addcovar <- as.matrix(addcovar)
  if(!is.null(intcovar)) intcovar <- as.matrix(intcovar)

  if(model=="2part") res <- matrix(ncol=3,nrow=n.perm)
  else res <- 1:n.perm

  if(verbose) { # if verbose, print out a tracing information
    # rnd: how often to print tracing information
    if(verbose > 1) rnd <- 1
    else {
      if(n.perm >= 1000) rnd <- 20
      else if(n.perm >= 100) rnd <- 5
      else rnd <- 1
    }
  }

  if(method=="mr-imp") # save version with missing genotypes 
    tempcross <- cross
  if(method=="mr-argmax") # impute genotypes
    cross <- fill.geno(cross,method="argmax")

  for(i in 1:n.perm) {
    if(verbose && i/rnd == round(i/rnd))
      cat("Permutation", i, "\n")

    # impute genotypes for method "mr-imp"
    if(method=="mr-imp") cross <- fill.geno(tempcross)

    o <- sample(1:n.ind)
    cross$pheno <- cross$pheno[o,,drop=FALSE]
    if(!is.null(addcovar)) addcovarp <- addcovar[o,,drop=FALSE]
    if(!is.null(intcovar)) intcovarp <- intcovar[o,,drop=FALSE]
    if(!is.null(weights)) weights <- weights[o]
    tem <- scanone(cross,,pheno.col,model,method,addcovarp,
                   intcovarp,weights,upper,ties.random,start,
                   maxit,tol,n.perm= -1)
    if(model=="2part")
      res[i,] <- apply(tem[,3:5], 2, max,na.rm=TRUE)
    else res[i] <- max(tem[,3],na.rm=TRUE)
  }

  if(model=="2part") colnames(res) <- c("LOD.p.mu", "LOD.p", "LOD.mu")

  attr(res,"method") <- method
  attr(res,"model") <- model
  attr(res,"type") <- class(cross)[1]
  res
}

##################################################################
# summarize one or more scanone results for given thresholds
##################################################################
summary.scanone <- function(..., threshold=0)
{
  # take the input objects
  object <- list(...)

  # if last thing in "..." is a number, take it as the threshold
  if(length(object) > 1 && missing(threshold) && is.numeric(object[[length(object)]])) {
    threshold <- object[[length(object)]]
    object <- object[-length(object)]
  }

  n.scanone <- length(object)
  
  ### data checking ... need to be more thorough
  # all scans need to be for the same chromosomes and markers, etc.
  if(n.scanone > 1) {
    chr <- object[[1]]$chr
    for(i in 2:n.scanone)
      if(length(object[[i]]$chr) != length(chr))
        stop("All scan results need to be for the same chromosomes and markers")
  }

  ### start to summarize
  # single scan result, this is easier
  if(n.scanone == 1) { # for a single scanone result, produce a longer summary
    object <- object[[1]]
    # pick off the maximum from each chromosome
    # this is made complicated to avoid returning multiple rows
    #     from any one chromosome
    out <- lapply(split(object,object[,1]),
                  function(b) {
                    d <- which(b[,3]==max(b[,3]))
                    if(length(d) > 1) d <- median(d)
                    b[d,] }) 
    results <- out[[1]]
    if(length(out) > 1)
      for(i in 2:length(out))
        results <- rbind(results,out[[i]])
    # deal with threshold ...
    # for one scanone result, threshold can be either a single value
    # or a vector corresponding to different significance levels
    # find the ones above threshold
    threshold <- sort(as.vector(threshold))
    n.sig <- length(threshold)
    # take the peaks
    idx <- which(results[,3] >= min(threshold))
    lod <- results[idx,3]
    results <- results[idx,]
    if(n.sig > 1) {
      # make significance code
      Sig <- rep("", length(idx))
      for(i in 1:length(idx)) {
        for(j in 1:n.sig) {
          if(lod[i]>threshold[j])
            Sig[i] <- paste(Sig[i],"*", sep="")
        }
      }
      results <- cbind(results, Sig)
    }
  }
  
  else { # multiple scanone results, report only chr, pos and lod
    #### deal with threshold ...
    # threshold should be a matrix where each row corresponds to a significance level
    # and each column corresponds to a scan result
    if(length(threshold) == 1) {
      # if threshold is a single value, make it a one row matrix 
      threshold <- matrix(rep(threshold, n.scanone), nrow=1)
    }
    else if(is.vector(threshold)) {
      # it's a vector. I assume it's for different scan results with one sig level
      if(length(threshold) != n.scanone)
        stop("The input threshold vector has wrong length")
      # convert to a one row matrix
      threshold <- matrix(threshold, nrow=1)
    }
    else if(is.matrix(threshold)) {
      # it's a matrix, rows are for sig levels and columns are for scan results
      # if it has one column, I assume the same threshold value will apply
      # to all scan results
      if(dim(threshold)[2] ==1)
        threshold <- matrix(rep(threshold, n.scanone), ncol=n.scanone)
      else if(dim(threshold)[2] != n.scanone)
        stop("The input threshold matrix has wrong dimension")
      # threshold should be in ascending order
      threshold <- apply(threshold, 2, sort)
    }
    else { # shouldn't come here
      stop("Wrong threshold")
    }
    
    # number of significance levels
    n.sig <- dim(threshold)[1]
    
    results <- NULL
    # chr and pos for the scan results
    chrpos <- object[[1]][c("chr","pos")]
    # put all lod scores together
    lod <- object[[1]]["lod"]
    for(i in 2:n.scanone)
      lod <- cbind(lod, object[[i]][,"lod"])
    # split chrpos and lod by chromosome
    chrpos.chr <- split(chrpos, chrpos[,"chr"])
    lod.chr <- split(lod,chrpos[,"chr"])
  
    # loop through the chromosomes - I could do this using lapply
    # but I'd rather use a loop for readability
    for(i in 1:length(lod.chr)) {
      ### find the peaks in this chromosome
      idx.peak <- NULL
      for(j in 1:n.scanone) {
        lodj <- lod.chr[[i]][,j]
        maxlodj <- max(lodj)
        # if this one didn't exceed any threshold, go to next one
        if(maxlodj < min(threshold[,j]))
          next;
        idx <- which(lodj==maxlodj)
        if(length(idx)>1) {
          # multiple peaks with the same lod,
          # pick one with the maximum total lod for all traits
          lod.tmp <- matrix(lod.chr[[i]][idx,-j], nrow=length(idx))
          lod.sum <- apply(lod.tmp, 1, sum)
          idx <- idx[which(lod.sum==max(lod.sum))]
        }
        idx.peak <- c(idx.peak, idx)
      }
      if(length(idx.peak) == 0) # no LOD exceed threshold on this chromosome, go to next one
        next
      else if(length(idx.peak) >= 1) { # if several scanone results peaked at the same position
        idx.peak <- unique(idx.peak)
      #  result <- rbind(result, lod.chr[[i]][idx.peak,])
      }
      # make summary information for this chromosome
      # -- this part of code is awkward because I need to add significance code...
      
      for(nn in 1:length(idx.peak)) {
        summary.tmp <- chrpos.chr[[i]][idx.peak[nn],]
        for(j in 1:n.scanone) {
          lod.tmp <- lod.chr[[i]][idx.peak[nn],j] # lod score for this scanresult on this peak
          # make sig code
          Sig <- ""
          for(k in 1:n.sig) {
            if(lod.tmp > threshold[k,j])
              Sig <- paste(Sig, "*", sep="")
          }
          summary.tmp <- cbind(summary.tmp, lod.tmp, Sig)
        }
        # make results
        results <- rbind(results, summary.tmp)
      }
    }
    # column names for summary result
    colnames(results) <- c("chr", "pos",
                          paste(rep(c("lod","Sig"), n.scanone), rep(1:n.scanone,each=2), sep=""))
  }

  # return 
  class(results) <- c("summary.scanone","data.frame")
  results
}

# print output of summary.scanone
print.summary.scanone <-
function(x,...)
{

  if(nrow(x) == 0) {
    cat("    There were no LOD peaks above the threshold.\n")
  }

  else {
    # find the column number of significance code
    idx.sig <- grep("Sig",colnames(x))
    # round the doubles
    x[,-c(1,2,idx.sig)] <- round(data.frame(x[,-c(1:2,idx.sig)]),6)
    cat("\n")
    print.data.frame(x,digits=2)
    cat("\n")
  }
}

# pull out maximum LOD peak, genome-wide
max.scanone <-
function(..., chr, na.rm=TRUE)
{
  dots <- list(...)[[1]]
  if(missing(chr)) {
    maxlod <- max(dots[,3],na.rm=TRUE)
    dots <- dots[!is.na(dots[,3]) & dots[,3]==maxlod,]
    return(summary.scanone(dots,threshold=0))
  }
  else {
    res <- NULL
    for(i in seq(along=chr)) {
      temp <- dots[dots[,1]==chr[i],]
      maxlod <- max(temp[,3],na.rm=TRUE)
      temp <- temp[!is.na(temp[,3]) & temp[,3]==maxlod,]
      res <- rbind(res,temp)
    }
    return(summary.scanone(res,threshold=0))
  }
}

# end of scanone.R
