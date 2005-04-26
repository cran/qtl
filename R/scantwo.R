######################################################################
#
# scantwo.R
#
# copyright (c) 2001-5, Karl W Broman, Johns Hopkins University,
#            Hao Wu, and Brian Yandell
# last modified Mar, 2005
# first written Nov, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Hao Wu (The Jackson Lab) wrote the initial code for the imputation
# method and summary.scantwo functions.  Brian Yandell made further
# modifications/enhancements to summary.scantwo, but Karl re-wrote
# most of it later.
#
# Part of the R/qtl package
# Contains: scantwo, scantwo.perm, summary.scantwo
#           print.summary.scantwo, max.scantwo
#
######################################################################

######################################################################
#
# scantwo: Do 2-dimensional genome scan with a two-QTL model,
#          calculating joint LOD scores and LOD scores testing
#          epistasis.
#
######################################################################

scantwo <-
function(cross, chr, pheno.col=1,
         method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
         model=c("normal","binary"),
         addcovar=NULL, intcovar=NULL, weights=NULL,
         run.scanone=TRUE, incl.markers=FALSE, maxit=4000, tol=1e-4,
         verbose=TRUE, n.perm)
{
  method <- match.arg(method)
  model <- match.arg(model)
  
  origcross <- cross

  # pull out chromosomes to be scanned
  if(!missing(chr)) cross <- subset(cross,chr=chr)
  if(missing(n.perm)) n.perm <- 0

  # check phenotypes and covariates; drop individuals with missing values
  # in case of permutation test, only do checks once
  if(n.perm>=0) { 
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
  chrtype <- sapply(cross$geno,class)

  if(model=="binary") {
    if(method != "em") {
      method <- "em"
      warning("Only EM algorithm coded for binary traits")
    }
    if(any(chrtype=="X")) {
      sexpgm <- getsex(cross)
      if(!is.null(sexpgm$sex) || !is.null(sexpgm$pgm)) {
        cross <- subset(cross,chr = chrtype!="X")
        n.chr <- nchr(cross)
        n.ind <- nind(cross)
        warning("X chromosome is not yet working for binary traits; dropping it.")
      }
    }
    if(!is.null(weights)) {
      weights <- NULL
      warning("weights ignored for binary traits.")
    }

    u <- unique(pheno)
    if(any(u!=0 & u!=1))
      stop("Phenotypes must be either 0 or 1.")
  }

  # Problems with EX w/ X chromosome: just use H-K for now.
  if(model=="normal" && any(chrtype=="X") && method=="em") {
    sexpgm <- getsex(cross)
    if(!is.null(sexpgm$sex) || !is.null(sexpgm$pgm)) {
      warning("EM is not yet working for X chromosomes; using H-K instead.")
      method <- "hk"
    }
  }

  # if n.perm specified, do a permutation test
  if(n.perm>0) { 
    return(scantwo.perm(cross, pheno.col, method, model, addcovar,
                        intcovar, weights, incl.markers,
                        maxit, tol, verbose, n.perm))
  }

  if(n.perm == 0) { # not in the midst of permutations
    if(method=="mr-argmax")
      cross <- fill.geno(cross,method="argmax")
    if(method=="mr-imp")
      cross <- fill.geno(cross,method="imp")
  }

  # weights of individuals
  if(model == "normal") {
    if(is.null(weights))
      weights <- rep(1, nind(cross))
    if(length(weights) != nind(cross))
      stop("weights should either be NULL or a vector of length n.ind")
    if(any(weights) <= 0)
      stop("weights should be entirely positive")
    weights <- sqrt(weights)
  }

  if(run.scanone) { # also do scanone
    if(verbose) cat(" --Running scanone\n")
    temp <- scanone(cross, pheno.col=pheno.col, method=method, model=model,
                    addcovar=addcovar, intcovar=intcovar, weights=weights,
                    maxit=maxit, tol=tol, verbose=FALSE)
    nam <- rownames(temp)
    out.scanone <- temp[,3]
    names(out.scanone) <- nam
    if(verbose) cat(" --Running scantwo\n")
  }

  if(method=="mr" || method=="mr-imp" || method=="mr-argmax") { # marker regression
    # number of genotypes on each chromosome, 
    #     combine the genetic maps for all chromosomes
    map <- unlist(pull.map(cross))
    names(map) <- unlist(lapply(pull.map(cross),names))
    n.pos <- nmar(cross)
    gmap <- data.frame(chr=rep(names(cross$geno),n.pos),
                       pos=map,
                       eq.spacing=rep(1,sum(n.pos)),
                       xchr=rep(sapply(cross$geno,class)=="X",nmar(cross)))

    # number of possible genotypes for each chromosome
    n.gen <- 1:n.chr
    for(i in 1:n.chr) { 
      if(chrtype[i]=="X") 
        sexpgm <- getsex(cross)
      else 
        sexpgm <- NULL

      gen.names <- getgenonames(type, chrtype[i], "full", sexpgm)
      n.gen[i] <- length(gen.names)
    }
  } # end of if(method=="mr")

  else { # all methods except "mr"
    # check for genotype probabilities or simulated genotypes
    steps <- rep(0,n.chr) # step length on each chromosome
    if(method=="imp") {
      for(i in 1:n.chr) {
        if(is.na(match("draws",names(cross$geno[[i]])))) {
          # need to run sim.geno
          warning("First running sim.geno.")
          cross <- sim.geno(cross)
        }
        steps[i] <- attr(cross$geno[[i]]$draws,"step")
      }

      # make sure all chromosomes have the same number of imputations
      n.draws <- sapply(cross$geno, function(a) dim(a$draws)[3])
      if(length(unique(n.draws)) > 1) {
        warning("Re-running sim.geno to have a fixed number of imputations.")
        cross <- sim.geno(cross, n.draws=max(n.draws),
                          step=attr(cross$geno[[1]]$draws,"step"),
                          off.end=attr(cross$geno[[1]]$draws,"off.end"))
      }
      n.draws <- max(n.draws)
    }
    else { # H-K or EM
      for(i in 1:n.chr) {
        if(is.na(match("prob",names(cross$geno[[i]])))) {
          # need to run calc.genoprob
          warning("First running calc.genoprob.")
          cross <- calc.genoprob(cross)
        }
        steps[i] <- attr(cross$geno[[i]]$prob,"step")
      }
    }

    # number of genotypes on each chromosome, 
    #     construct the genetic map for all chromosomes
    #     and possibly drop marker positions
    gmap <- NULL
    n.pos <- n.gen <- rep(0,n.chr) 
    keep.pos <- vector("list",n.chr)
    some.dropped <- rep(FALSE,n.chr)

    for(i in 1:n.chr) { 
      if(chrtype[i]=="X") 
        sexpgm <- getsex(cross)
      else 
        sexpgm <- NULL

      gen.names <- getgenonames(type, chrtype[i], "full", sexpgm)
      n.gen[i] <- length(gen.names)

      # construct the genetic map for this chromesome
      if(method=="imp") 
        map <- create.map(cross$geno[[i]]$map,
                          attr(cross$geno[[i]]$draws,"step"),
                          attr(cross$geno[[i]]$draws,"off.end"))
      else
        map <- create.map(cross$geno[[i]]$map,
                          attr(cross$geno[[i]]$prob,"step"),
                          attr(cross$geno[[i]]$prob,"off.end"))

      if(is.matrix(map)) map <- map[1,] # in case of sex-specific map
  
      w <- names(map)
      o <- grep("^loc\-*[0-9]+",w)

      if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
        w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
      map <- cbind(chr=rep(names(cross$geno)[i],length(map)),
                   pos=as.data.frame(map) )
      rownames(map) <- w 

      # equally spaced positions
      if(steps[i]==0)  # just use markers
        eq.sp.pos <- rep(1,nrow(map))
      else {
        eq.sp.pos <- seq(min(map[,2]),max(map[,2]),by=steps[i])
        wh.eq.sp <- match(eq.sp.pos,map[,2])
        if(any(is.na(wh.eq.sp))) { # this shouldn't happen
          warning("Possible error in determining the equally spaced positions.")
          wh.eq.sp <- wh.eq.sp[!is.na(wh.eq.sp)]
        }
        eq.sp.pos <- rep(0,nrow(map))
        eq.sp.pos[wh.eq.sp] <- 1
      }
      if(!incl.markers && any(eq.sp.pos==0)) {
        keep.pos[[i]] <- (seq(along=eq.sp.pos))[eq.sp.pos==1]
        map <- map[eq.sp.pos==1,]
        eq.sp.pos <- eq.sp.pos[eq.sp.pos==1]
        some.dropped[i] <- TRUE # indicates some positions were dropped
      }
      else keep.pos[[i]] <- seq(along=eq.sp.pos)
      gmap <- rbind(gmap, cbind(map,eq.spacing=eq.sp.pos,
                                xchr=(class(cross$geno[[i]])=="X")))
      n.pos[i] <- length(keep.pos[[i]])

      # Revise X chromosome genotype probabilities or imputations 
      if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss")) {
        if(method=="imp") 
          cross$geno[[i]]$draws <-
            reviseXdata(type, "full", sexpgm, draws=cross$geno[[i]]$draws)
        else if(method=="hk" || method=="em") {
          oldXchr <- subset(cross, chr=i)
          cross$geno[[i]]$prob <-
            reviseXdata(type, "full", sexpgm, prob=cross$geno[[i]]$prob)
        }
        else 
          cross$geno[[i]]$data <-
            reviseXdata(type, "full", sexpgm, data=cross$geno[[i]]$data)
      }

    } # end loop over chromosomes
  } # end of if/else for method="mr" vs other 

  # columns in result matrix for each chromosome
  wh.col <- vector("list",n.chr)
  first.pos <- cumsum(c(1,n.pos))
  for(i in 1:n.chr)
    wh.col[[i]] <- seq(first.pos[i],by=1,length=n.pos[i])

  # initialize the results matrix
  results <- matrix(0,ncol=sum(n.pos), nrow=sum(n.pos))

  # do the 2-dimensional genome scan
  for(i in 1:n.chr) { # loop over the 1st chromosome
    for(j in i:n.chr) { # loop over the 2nd chromosome

      if(chrtype[i]=="X" || chrtype[j]=="X") {
        ac <- revisecovar(sexpgm,addcovar)
        n.ac <- ifelse(is.null(ac),0,ncol(ac))
        ic <- revisecovar(sexpgm,intcovar)
        n.ic <- ifelse(is.null(ic),0,ncol(ic))
      }
      else {
        ac <- addcovar
        n.ac <- n.addcovar
        ic <- intcovar
        n.ic <- n.intcovar
      }

      # print the current working pair
      if(verbose) cat(paste(" (", names(cross$geno)[i], ",",
                          names(cross$geno)[j],")\n",sep=""))

      if(method=="imp") {
        z <- .C("R_scantwo_imp",
                as.integer(n.ind),
                as.integer(i==j),
                as.integer(n.pos[i]),
                as.integer(n.pos[j]),
                as.integer(n.gen[i]),
                as.integer(n.gen[j]),
                as.integer(n.draws),
                as.integer(cross$geno[[i]]$draws[,keep.pos[[i]],]),
                as.integer(cross$geno[[j]]$draws[,keep.pos[[j]],]),
                as.double(ac),
                as.integer(n.ac),
                as.double(ic),
                as.integer(n.ic),
                as.double(pheno),
                as.double(weights),
                result=as.double(rep(0,2*n.pos[i]*n.pos[j])),
                PACKAGE="qtl")
        z <- array(z$result,dim=c(n.pos[i], n.pos[j], 2)) # rearrange the result 

        # update the final result matrix
        results[wh.col[[i]],wh.col[[j]]] <- z[,,1]
        if(i != j) results[wh.col[[j]],wh.col[[i]]] <- t(z[,,2])
        else { # do this just once: do null model and get neg log10 likelihood
          if(i==1) { 
            if(n.ac > 0)
              resid0 <- lm(pheno ~ ac, weights=weights^2)$resid
            else
              resid0 <- lm(pheno ~ 1, weights=weights^2)$resid
            sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
            nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
          }
        }
      }
      else if(model=="normal" && (method=="hk" || method=="em")) {
        if(i==j) { # same chromosome

          if(i==1) { # first time! do null model and get neg log10 likelihood
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


          if(verbose>1) cat("  --Calculating joint probs.\n")

          if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss")) {
            # calculate joint genotype probabilities for all pairs of positions
            stp <- attr(oldXchr$geno[[1]]$prob, "step")
            oe <- attr(oldXchr$geno[[1]]$prob, "off.end")
            err <- attr(oldXchr$geno[[1]]$prob, "error.prob")
            mf <- attr(oldXchr$geno[[1]]$prob, "map.function")

            temp <- calc.pairprob(oldXchr,stp,oe,err,mf)
          }
          else {
            # calculate joint genotype probabilities for all pairs of positions
            stp <- attr(cross$geno[[i]]$prob, "step")
            oe <- attr(cross$geno[[i]]$prob, "off.end")
            err <- attr(cross$geno[[i]]$prob, "error.prob")
            mf <- attr(cross$geno[[i]]$prob, "map.function")

            temp <- calc.pairprob(subset(cross,chr=i),stp,oe,err,mf)
          }

          # pull out positions from genotype probs
          if(some.dropped[i]) {
            # figure out pos'ns corresponding to columns of temp
            nc <- ncol(cross$geno[[i]]$prob)
            ind <- matrix(rep(1:nc,nc),ncol=nc)
            w <- lower.tri(ind)
            ind <- cbind(first=t(ind)[w],second=ind[w])

            # which part to keep
            keep <- apply(ind,1,function(a,b) all(!is.na(match(a,b))),
                          keep.pos[[i]])
            temp <- temp[,keep,,]
          }

          # revise pair probilities for X chromosome
          if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss")) {
            temp <- reviseXdata(type, "full", sexpgm, pairprob=temp)
            temp[temp < 1e-5] <- 1e-5 # << temp fix for problems with X chromosome
          }

          if(verbose>1) cat("  --Done.\n")

          if(method=="hk") 
            z <- .C("R_scantwo_1chr_hk", 
                    as.integer(n.ind),
                    as.integer(n.pos[i]),
                    as.integer(n.gen[i]),
                    as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                    as.double(temp),
                    as.double(ac),
                    as.integer(n.ac),
                    as.double(ic),
                    as.integer(n.ic),
                    as.double(pheno),
                    as.double(weights),
                    result=as.double(rep(0,n.pos[i]^2)),
                    PACKAGE="qtl")
          else
            z <- .C("R_scantwo_1chr_em", 
                    as.integer(n.ind),
                    as.integer(n.pos[i]),
                    as.integer(n.gen[i]),
                    as.double(temp),
                    as.double(ac),
                    as.integer(n.ac),
                    as.double(ic),
                    as.integer(n.ic),
                    as.double(pheno),
                    as.double(weights),
                    result=as.double(rep(0,n.pos[i]^2)),
                    as.integer(maxit),
                    as.double(tol),
                    as.integer(verbose),
                    PACKAGE="qtl")

          rm(temp) # remove the joint genotype probabilities

          # re-organize results
          results[wh.col[[i]],wh.col[[i]]] <-
            matrix(z$result,ncol=n.pos[i])
        } # end same chromosome
        else {
          if(method=="hk")
            z <- .C("R_scantwo_2chr_hk",
                    as.integer(n.ind),
                    as.integer(n.pos[i]),
                    as.integer(n.pos[j]),
                    as.integer(n.gen[i]),
                    as.integer(n.gen[j]),
                    as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                    as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                    as.double(ac),
                    as.integer(n.ac),
                    as.double(ic),
                    as.integer(n.ic),
                    as.double(pheno),
                    as.double(weights),
                    full=as.double(rep(0,n.pos[i]*n.pos[j])),
                    int=as.double(rep(0,n.pos[i]*n.pos[j])),
                    PACKAGE="qtl")
          else 
            z <- .C("R_scantwo_2chr_em",
                    as.integer(n.ind),
                    as.integer(n.pos[i]),
                    as.integer(n.pos[j]),
                    as.integer(n.gen[i]),
                    as.integer(n.gen[j]),
                    as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                    as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                    as.double(ac),
                    as.integer(n.ac),
                    as.double(ic),
                    as.integer(n.ic),
                    as.double(pheno),
                    as.double(weights),
                    full=as.double(rep(0,n.pos[i]*n.pos[j])),
                    int=as.double(rep(0,n.pos[i]*n.pos[j])),
                    as.integer(maxit),
                    as.double(tol),
                    as.integer(verbose),
                    PACKAGE="qtl")

          results[wh.col[[j]],wh.col[[i]]] <-
            t(matrix(z$full,ncol=n.pos[j]))
          results[wh.col[[i]],wh.col[[j]]] <-
            matrix(z$int,ncol=n.pos[j])
        } # end same chromosome
      }
      else if(model=="binary" && method=="em") {
        if(i==j) { # same chromosome

          if(i==1) { # first time! do null model and get neg log10 likelihood
            if(n.ac > 0)
              nullfit <- glm(pheno ~ ac, family=binomial(link=logit))
            else
              nullfit <- glm(pheno ~ 1, family=binomial(link=logit))
            fitted <- nullfit$fitted
            nullcoef <- nullfit$coef
            nllik0 <- -sum(pheno*log10(fitted) + (1-pheno)*log10(1-fitted))
            if(verbose > 1) cat("null log lik: ", nllik0, "\n")
          }

          start <- c(rep(nullcoef[1],n.gen[i]),rep(0,n.gen[i]-1),
                     nullcoef[-1],rep(0,n.gen[i]*n.gen[i]+
                                      (n.gen[i]*(n.gen[i]-1)*n.ic)))

          if(verbose>1) cat("  --Calculating joint probs.\n")

          if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss")) {
            # calculate joint genotype probabilities for all pairs of positions
            stp <- attr(oldXchr$geno[[1]]$prob, "step")
            oe <- attr(oldXchr$geno[[1]]$prob, "off.end")
            err <- attr(oldXchr$geno[[1]]$prob, "error.prob")
            mf <- attr(oldXchr$geno[[1]]$prob, "map.function")

            temp <- calc.pairprob(oldXchr,stp,oe,err,mf)
          }
          else {
            # calculate joint genotype probabilities for all pairs of positions
            stp <- attr(cross$geno[[i]]$prob, "step")
            oe <- attr(cross$geno[[i]]$prob, "off.end")
            err <- attr(cross$geno[[i]]$prob, "error.prob")
            mf <- attr(cross$geno[[i]]$prob, "map.function")

            temp <- calc.pairprob(subset(cross,chr=i),stp,oe,err,mf)
          }

          # pull out positions from genotype probs
          if(some.dropped[i]) {
            # figure out pos'ns corresponding to columns of temp
            nc <- ncol(cross$geno[[i]]$prob)
            ind <- matrix(rep(1:nc,nc),ncol=nc)
            w <- lower.tri(ind)
            ind <- cbind(first=t(ind)[w],second=ind[w])

            # which part to keep
            keep <- apply(ind,1,function(a,b) all(!is.na(match(a,b))),
                          keep.pos[[i]])
            temp <- temp[,keep,,]
          }

          # revise pair probilities for X chromosome
          if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss")) {
            temp <- reviseXdata(type, "full", sexpgm, pairprob=temp)
            temp[temp < 1e-5] <- 1e-5 # << temp fix for problems with X chromosome
          }

          if(verbose>1) cat("  --Done.\n")

          z <- .C("R_scantwo_1chr_binary_em", 
                  as.integer(n.ind),
                  as.integer(n.pos[i]),
                  as.integer(n.gen[i]),
                  as.double(temp),
                  as.double(ac),
                  as.integer(n.ac),
                  as.double(ic),
                  as.integer(n.ic),
                  as.integer(pheno),
                  as.double(start),
                  result=as.double(rep(0,n.pos[i]^2)),
                  as.integer(maxit),
                  as.double(tol),
                  as.integer(verbose),
                  PACKAGE="qtl")

          rm(temp) # remove the joint genotype probabilities

          # re-organize results
          results[wh.col[[i]],wh.col[[i]]] <-
            matrix(z$result,ncol=n.pos[i])
        } # end same chromosome
        else {
          start <- c(rep(nullcoef[1],n.gen[i]),rep(0,n.gen[j]-1),
                     nullcoef[-1],rep(0,n.gen[i]*n.gen[j]+
                                      (n.gen[i]*(n.gen[j]-1)*n.intcovar)));

          z <- .C("R_scantwo_2chr_binary_em",
                  as.integer(n.ind),
                  as.integer(n.pos[i]),
                  as.integer(n.pos[j]),
                  as.integer(n.gen[i]),
                  as.integer(n.gen[j]),
                  as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                  as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                  as.double(ac),
                  as.integer(n.ac),
                  as.double(ic),
                  as.integer(n.ic),
                  as.integer(pheno),
                  as.double(start),
                  full=as.double(rep(0,n.pos[i]*n.pos[j])),
                  int=as.double(rep(0,n.pos[i]*n.pos[j])),
                  as.integer(maxit),
                  as.double(tol),
                  as.integer(verbose),
                  PACKAGE="qtl")

          results[wh.col[[j]],wh.col[[i]]] <-
            t(matrix(z$full,ncol=n.pos[j]))
          results[wh.col[[i]],wh.col[[j]]] <-
            matrix(z$int,ncol=n.pos[j])
        } # end same chromosome
      }
      else { # marker regression
        # replace missing and partially informative genotypes with 0's
        datai <- cross$geno[[i]]$data
        datai[is.na(datai)] <- 0
        if(type=="f2" || type=="f2ss") datai[datai>3] <- 0
        else if(type=="4way") datai[datai>4] <- 0

        if(i==j) { # same chromosome

          z <- .C("R_scantwo_1chr_mr",
                  as.integer(n.ind),
                  as.integer(n.pos[i]),
                  as.integer(n.gen[i]),
                  as.integer(datai),
                  as.double(ac),
                  as.integer(n.ac),
                  as.double(ic),
                  as.integer(n.ic),
                  as.double(pheno),
                  as.double(weights),
                  result=as.double(rep(0,n.pos[i]^2)),
                  PACKAGE="qtl")

          # re-organize results
          results[wh.col[[i]],wh.col[[i]]] <-
            matrix(z$result,ncol=n.pos[i])
        } # end same chromosome
        else {
          
          # replace missing and partially informative genotypes with 0's
          dataj <- cross$geno[[j]]$data
          dataj[is.na(dataj)] <- 0
          if(type=="f2" || type=="f2ss") dataj[dataj>3] <- 0
          else if(type=="4way") dataj[dataj>4] <- 0

          z <- .C("R_scantwo_2chr_mr",
                  as.integer(n.ind),
                  as.integer(n.pos[i]),
                  as.integer(n.pos[j]),
                  as.integer(n.gen[i]),
                  as.integer(n.gen[j]),
                  as.integer(datai),
                  as.integer(dataj),
                  as.double(ac),
                  as.integer(n.ac),
                  as.double(ic),
                  as.integer(n.ic),
                  as.double(pheno),
                  as.double(weights),
                  full=as.double(rep(0,n.pos[i]*n.pos[j])),
                  int=as.double(rep(0,n.pos[i]*n.pos[j])),
                  PACKAGE="qtl")

          results[wh.col[[j]],wh.col[[i]]] <-
            t(matrix(z$full,ncol=n.pos[j]))
          results[wh.col[[i]],wh.col[[j]]] <-
            matrix(z$int,ncol=n.pos[j])
        } # end same chromosome
      }
    
    } # end loop over second chr
  } # end loop over first chromosome

  if(method=="hk" || method=="em") # subtr null neg log lik from lower tri
    results[lower.tri(results)] <- nllik0 - results[lower.tri(results)]

  # If the X chromosome was included, need to do an adjustment...
  scanoneX <- NULL
  if(any(gmap[,4])) { # the X chromosome was included

    # determine which covariates belong in null hypothesis
    temp <- scanoneXnull(type, sexpgm)
    adjustX <- temp$adjustX
    dfX <- temp$dfX
    sexpgmcovar <- temp$sexpgmcovar
      
    if(adjustX) {
      if(method=="mr") {
        warning("The X chr may not be working properly for scantwo with method mr.") 
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
          sigX <- sqrt(sum((residX*weights)^2)/n.ind)
          nllikX <- -sum(dnorm(residX,0,sigX/weights,log=TRUE))/log(10)
        }

        wh <- ((gmap[row(results),4] | gmap[col(results),4]) & lower.tri(results))
        results[wh] <- results[wh] + nllikX - nllik0

        if(run.scanone) {
          if(verbose) cat(" --Running scanone with special X chr covariates\n")

          notxchr <- which(sapply(cross$geno,class)!="X")
          if(length(notxchr) > 0) {
            temp <- scanone(subset(cross,chr=notxchr),
                            pheno.col=pheno.col, method=method,
                            addcovar=cbind(ac,sexpgmcovar),
                            intcovar=ic, weights=weights,
                            maxit=maxit, tol=tol, verbose=FALSE)


            nam <- rownames(temp)
            scanoneX <- temp[,3]
            names(scanoneX) <- nam

            scanoneX <- c(scanoneX,out.scanone[rownames(gmap)][gmap[,4]])
            scanoneX <- scanoneX[rownames(gmap)]
          }
          else {
            scanoneX <- out.scanone[rownames(gmap)][gmap[,4]]
            scanoneX <- scanoneX[rownames(gmap)]
          }
        }
      }
    }
  }
  

  if(any(is.na(results) | results < -1e-6 | results == Inf))
    warning("Some LOD scores NA, Inf or < 0")
  
  # output has 2 fields, lod and map
  out <- list(lod=results,map=gmap,scanoneX=scanoneX)
  class(out) <- "scantwo"

  if(run.scanone) # also did scanone
    diag(out$lod) <- out.scanone[rownames(out$map)]

  attr(out,"method") <- method
  attr(out,"type") <- type
  out
}

######################################################################
#
# scantwo.perm: Permutation test of scantwo
#
######################################################################

scantwo.perm <-
function(cross, pheno.col=1,
         method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
         model=c("normal","binary"),
         addcovar=NULL, intcovar=NULL, weights=NULL,
         incl.markers=FALSE, maxit=4000, tol=1e-4, verbose=FALSE,
         n.perm=1000) 
{
  method <- match.arg(method)
  model <- match.arg(model)

  n.ind <- nind(cross)
  addcovarp <- intcovarp <- NULL
  if(!is.null(addcovar)) addcovar <- as.matrix(addcovar)
  if(!is.null(intcovar)) intcovar <- as.matrix(intcovar)

  if(method=="mr-imp") # save version with missing genotypes 
    tempcross <- cross
  if(method=="mr-argmax") # impute genotypes
    cross <- fill.geno(cross,method="argmax")

  # initialize the result matrix
  # the first row is for full model comparison
  # the second row is for additive model comparison
  res <- matrix(ncol=2,nrow=n.perm)
  for(i in 1:n.perm) {
    if(verbose) cat("Permutation", i, "\n")

    # impute genotypes for method "mr-imp"
    if(method=="mr-imp") cross <- fill.geno(tempcross)

    o <- sample(1:n.ind)
    cross$pheno <- cross$pheno[o,,drop=FALSE]
    if(!is.null(addcovar)) addcovarp <- addcovar[o,,drop=FALSE]
    if(!is.null(intcovar)) intcovarp <- intcovar[o,,drop=FALSE]
    tem <- scantwo(cross,  pheno.col=pheno.col,
                   method=method, model=model, addcovar=addcovarp,
                   intcovar=intcovarp, incl.markers=incl.markers,
                   weights=weights, run.scanone=FALSE, maxit=maxit,
                   tol=tol, verbose=FALSE, n.perm = -1)

    # take max of the two triangles
    res[i,1] <- max( tem$lod[tem$lod < Inf & row(tem$lod)>col(tem$lod)], na.rm=TRUE )
    res[i,2] <- max( tem$lod[tem$lod < Inf & row(tem$lod)<col(tem$lod)], na.rm=TRUE )
  }
  colnames(res) <- c("LOD.jnt","LOD.interxn")
  attr(res,"method") <- method
  res
}


######################################################################
#
# summerize the result from scantwo
#
######################################################################
summary.scantwo <-
function (object, thresholds = c(0, 0, 0),
          type = c("joint","interaction"), ...)
{
  type <- match.arg(type)
  
  if(length(thresholds) < 3) {
    if(length(thresholds) == 1) thresholds <- c(thresholds, 0, 0)
    else stop("You must give three thresholds: full, interaction and main\n")
  }

  thrfull <- thresholds[1]
  thrint <- thresholds[2]
  thrcond <- thresholds[3]

  lod <- object$lod
  map <- object$map

  # backward compatibility for previous version of R/qtl
  if(is.na(match("scanoneX",names(object)))) {
    warning("It would be best to re-run scantwo() with the R/qtl version 0.98 or later.")
    scanoneX <- NULL
  }
  else scanoneX <- object$scanoneX

  # deal with bad LOD score values
  if(any(is.na(lod) | lod < -1e-06 | lod == Inf)) 
    warning("Some LOD scores NA, Inf or < 0; set to 0")
  lod[is.na(lod) | lod < 0 | lod == Inf] <- 0

  # if there's no mainscan result, ignore the thresholds
  #     and don't include the 4 conditional LOD columns
  if(all(is.na(diag(lod)) | diag(lod) < 1e-10)) 
    includes.scanone <- FALSE
  else includes.scanone <- TRUE

  # If scanone results available, calculate conditional LOD scores
  if(includes.scanone) {
    d <- diag(lod)
    q1 <- matrix(rep(d,length(d)),ncol=length(d))
    q2 <- matrix(rep(d,length(d)),ncol=length(d),byrow=TRUE)

    if(!is.null(scanoneX) && any(map[,4])) {
      d <- scanoneX
      q1X <- matrix(rep(d,length(d)),ncol=length(d))
      q2X <- matrix(rep(d,length(d)),ncol=length(d),byrow=TRUE)
      q1[map[,4],] <- q1X[map[,4],]
      q2[,map[,4]] <- q2X[,map[,4]]
    }

    q1[lower.tri(q1)] <- t(q2)[lower.tri(q2)]
    condlod <- abs(lod - t(lod)) - q1
    diag(condlod) <- 0
  }
  else condlod <- NULL

  # Negative thresholds are interpreted relative to the maximum LOD score
  if(thrfull < 0) 
    thrfull <- max(0,max(lod[lower.tri(lod)]) + thrfull)
  if(thrint < 0) 
    thrint <- max(0,max(lod[upper.tri(lod)]) + thrint)
  if(thrcond < 0 && includes.scanone)
    thrcond <- max(0,max(condlod) + thrcond)
  
  crosstype <- attr(object, "type")
  if(is.null(crosstype)) {
    warning("No type attribute in input data; assuming backcross.")
    crosstype <- "bc"
  }

  # calculate the degree of freedom
  if(crosstype == "bc" || crosstype == "riself" || crosstype == 
      "risib") {
    df.int <- 1
    df.add <- 1
  }
  else if(crosstype == "f2") {
    df.int <- 4
    df.add <- 2
  }
  else if(crosstype == "4way") {
    df.int <- 9
    df.add <- 3
  }
  else {
    stop("Don't know what to do with cross type ", crosstype)
  }

  # chromsomes in the result
  chr <- unique(map[, 1])
  n.chr <- length(chr)

  # calculate the locations of each chromosome within the LOD matrix
  wh.index <- vector("list", n.chr)
  n <- nrow(map)
  for(i in 1:n.chr)
    wh.index[[i]] <- which(map[, 1] == chr[i])

  results <- NULL

  # go through each pair of chromosomes
  for(i in 1:n.chr) {
    for(j in i:n.chr) { 
      tmplod1 <- lod[wh.index[[j]], wh.index[[i]]]
      if(!is.null(condlod)) {
        if(i==j) tmpcondlod <- condlod[wh.index[[i]],wh.index[[i]]]
        else {
          tmpcondlod1 <- condlod[wh.index[[j]],wh.index[[i]]]
          tmpcondlod2 <- condlod[wh.index[[i]],wh.index[[j]]]
        }
      }

      if(i != j) tmplod2 <- lod[wh.index[[i]], wh.index[[j]]]
      else tmplod2 <- tmplod1
      

      if(type == "joint") {
        if(i == j) {
          tri <- lower.tri(tmplod1)
          lod.joint <- max(tmplod1[tri])
          idx <- which(tmplod1 == lod.joint & tri, arr.ind=TRUE)
        }
        else {
          lod.joint <- max(tmplod1)
          idx <- which(tmplod1 == lod.joint, arr.ind=TRUE)
        }
        if(nrow(idx)>1) idx <- idx[sample(nrow(idx),1),]
        idx.row <- idx[1]
        idx.col <- idx[2]
        
        lod.int <- tmplod2[idx.col, idx.row]
      }
      else { # interaction lod
        if(i == j) {
          tri <- upper.tri(tmplod2)
          lod.int <- max(tmplod2[tri])
          idx <- which(tmplod2 == lod.int & tri, arr.ind=TRUE)
        }
        else {
          lod.int <- max(tmplod2)
          idx <- which(tmplod2 == lod.int)
        }
        if(nrow(idx)>1) idx <- idx[sample(nrow(idx),1),]
        idx.row <- idx[2]
        idx.col <- idx[1]
        
        lod.joint <- tmplod1[idx.row, idx.col]
      }
      
      full.idx.row <- idx.row + wh.index[[j]][1] - 1
      full.idx.col <- idx.col + wh.index[[i]][1] - 1

      flag <- FALSE # a flag to indicate whether there's any peak on this pair
      if(lod.joint >= thrfull) {
        if(includes.scanone) {
          if(i==j) {
            lod.q1 <- tmpcondlod[idx.row,idx.col]
            lod.q2 <- tmpcondlod[idx.col,idx.row]
          }
          else {
            lod.q1 <- tmpcondlod1[idx.row,idx.col]
            lod.q2 <- tmpcondlod2[idx.col,idx.row]
          }
          
          if(lod.int >= thrint || min(c(lod.q1, lod.q2)) >= thrcond) {
            flag <- TRUE
            i.pos <- map[full.idx.col, 2]
            j.pos <- map[full.idx.row, 2]
            results <- rbind(results,
                             data.frame(chr[i], chr[j], i.pos, j.pos,
                                        lod.joint, 1 - pchisq(2 * log(10) * lod.joint,
                                                              df.int + 2 * df.add), 
                                        lod.int, 1 - pchisq(2 * log(10) * lod.int, df.int),
                                        lod.q1, 1 - pchisq(2 * log(10) * lod.q1, df.add),
                                        lod.q2, 1 - pchisq(2 * log(10) * lod.q2, df.add))
                             )
          }
        } 
        else { # no scanone output
          flag <- TRUE
          i.pos <- map[full.idx.col, 2]
          j.pos <- map[full.idx.row, 2]
          results <- rbind(results,
                           data.frame(chr[i], chr[j], i.pos, j.pos,
                                      lod.joint, 1 - pchisq(2 * log(10) * lod.joint,
                                                            df.int + 2 * df.add), 
                                      lod.int, 1 - pchisq(2 * log(10) * lod.int, df.int))
                           )
        }
        # give the new row (if any) a name
        if(flag) {
          mname <- rownames(map)
          rownames(results)[nrow(results)] <- paste(mname[full.idx.col], ":",
                                                    mname[full.idx.row], sep="")
        }
      } # lod joint above threshold

    } # end loop over chromosomes
  }

  if(is.null(results)) {
    results <- numeric(0)
  }
  else {
    if(includes.scanone) 
      colnames(results) <- c("chr1", "chr2", "pos1", "pos2", 
                             "lod.joint", "p.joint", "lod.int", "p.int", "lod.q1", 
                             "p.q1", "lod.q2", "p.q2")
    else colnames(results) <- c("chr1", "chr2", "pos1", "pos2", 
                                "lod.joint", "p.joint", "lod.int", "p.int")
    results <- as.data.frame(results)
  }
  class(results) <- c("summary.scantwo", "data.frame")
  results
}


print.summary.scantwo <-
function(x,...)
{
  if(length(x)==0) {
    cat("    There were no pairs of loci meeting the criteria.\n")
    invisible(return(NULL))
  }

  # column names
  cnames <- c("pos1", "pos2", "  LODjnt", "-logP",
              "  LODint", "-logP", "  LODq1", "-logP",
              "  LODq2", "-logP")

  # chr names
  chr1 <- paste("c",x[,1],sep="")
  chr2 <- paste("c",x[,2],sep="")

  # pad chr names with spaces; this isn't really necessary
  nchar.c1 <- nchar(chr1); max.nchar.c1 <- max(nchar.c1)
  nchar.c2 <- nchar(chr2); max.nchar.c2 <- max(nchar.c2)
  if(any(nchar.c1 < max.nchar.c1 | nchar.c2 < max.nchar.c2)) {
    for(i in 1:length(nchar.c2)) {
      if(nchar.c1[i] < max.nchar.c1)
        chr1[i] <- paste(paste(rep(" ", max.nchar.c1-nchar.c1[i]),collapse=""),
                         chr1[i],sep="")
      if(nchar.c2[i] < max.nchar.c2)
        chr2[i] <- paste(paste(rep(" ", max.nchar.c2-nchar.c2[i]),collapse=""),
                         chr2[i],sep="")
    }
  }
  chr <- paste(chr1,chr2,sep=":")

  # round the rest; take -log10(P-values)
  for(j in 3:ncol(x)) {
    if(j<5)
      x[,j] <- round(x[,j])
    else if(j %% 2)  # odd
      x[,j] <- round(x[,j],2)
    else
      x[,j] <- -round(log10(x[,j]),1)
  }

  res <- as.data.frame(x[,-(1:2)])
  names(res) <- cnames[1:ncol(res)]
  rownames(res) <- chr

  cat("\n")
  print.data.frame(res)
  cat("\n")
}

######################################################################
#
# max.scantwo:  Give maximum joint and intxnLODs for results of the
#               scantwo function
#
######################################################################

max.scantwo <-
function(..., na.rm=TRUE)
{
  dots <- list(...)[[1]]
  lod <- dots$lod
  map <- dots$map

  lod[is.na(lod) | lod == Inf | lod == -Inf] <- 0

  # maximum LODs
  max.jnt <- max(lod[row(lod)>col(lod)],na.rm=na.rm)
  max.int <- max(lod[row(lod)<col(lod)],na.rm=na.rm)

  # "zero" out everything but the maxima
  minmax <- c(min(max.jnt,max.int)/2)
  lod[row(lod)>col(lod) & !is.na(lod) &
      (lod<max.jnt & t(lod)<max.int)] <- minmax/10
  lod[row(lod)<col(lod) & !is.na(lod) &
      (t(lod)<max.jnt & lod<max.int)] <- minmax/10
  diag(lod) <- 0
  dots$lod <- lod

  # get locations of just the maxima
  summary(dots, c(minmax,0,0))
}


# end of scantwo.R
