######################################################################
#
# scantwo.R
#
# copyright (c) 2001-6, Karl W Broman, Johns Hopkins University,
#            Hao Wu, and Brian Yandell
#
# last modified Jun, 2006
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
#           scantwo.perm.engine
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
         model=c("normal","binary"),
         method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
         addcovar=NULL, intcovar=NULL, weights=NULL,
         run.scanone=TRUE, incl.markers=FALSE, maxit=4000, tol=1e-4,
         verbose=TRUE, n.perm, n.node=1)
{
  method <- match.arg(method)
  model <- match.arg(model)
  
  origcross <- cross
  
  # pull out chromosomes to be scanned
  if(!missing(chr)) cross <- subset(cross,chr=chr)
  if(missing(n.perm)) n.perm <- 0

  if(any(pheno.col < 1 | pheno.col > nphe(cross)))
    stop("pheno.col values should be between 1 and the no. phenotypes")

  # multiple phenotype for methods other than imp and hk
  if(length(pheno.col)>1 && n.perm <= 0 &&
     method!="imp" && method != "hk" ) {
    n.phe <- length(pheno.col)
    if(verbose) cat(" -Phenotype 1\n")
    output <- scantwo(cross, , pheno.col[1], model, method, addcovar,
                      intcovar, weights, run.scanone, incl.markers, maxit,
                      tol, verbose, n.perm, n.node)
    temp <- array(dim=c(nrow(output$lod), ncol(output$lod), n.phe))
    temp[,,1] <- output$lod
    output$lod <- temp
    for(i in 2:n.phe) {
      if(verbose) cat(" -Phenotype ", i, "\n")
      temp <- scantwo(cross, , pheno.col[i], model, method, addcovar,
                      intcovar, weights, run.scanone, incl.markers, maxit,
                      tol, verbose, n.perm, n.node)
      output$lod[,,i] <- temp$lod
      if(!is.null(output$scanoneX))
        output$scanoneX <- cbind(output$scanoneX, temp$scanoneX)
    }
    return(output)
  }


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
  else {  # come here only in permutation tests
    pheno <- as.matrix(cross$pheno[,pheno.col])
    if(is.null(addcovar)) n.addcovar <- 0
    else n.addcovar <- ncol(addcovar)
    if(is.null(intcovar)) n.intcovar <- 0
    else n.intcovar <- ncol(intcovar)
  }
  n.chr <- nchr(cross)
  n.ind <- nind(cross)
  n.phe <- length(pheno.col)
  type <- class(cross)[1]
  chrtype <- sapply(cross$geno,class)

  if(model=="binary") {
    if(method != "em") {
      method <- "em"
      warning("Only EM algorithm coded for binary traits")
    }
    if(any(chrtype=="X")) {
      sexpgm <- getsex(cross)
#      if(!is.null(sexpgm$sex) || !is.null(sexpgm$pgm)) {
#        cross <- subset(cross,chr = chrtype!="X")
#        n.chr <- nchr(cross)
#        n.ind <- nind(cross)
#        warning("X chromosome is not yet working for binary traits; dropping it.")
#      }
    }

    if(!is.null(weights)) {
      weights <- NULL
      warning("weights ignored for binary traits.")
    }

    u <- unique(pheno)
    if(any(u!=0 & u!=1))
      stop("Phenotypes must be either 0 or 1.")
  }

  # if n.perm specified, do a permutation test
  if(n.perm>0) { 
    return(scantwo.perm(cross, pheno.col, model, method, addcovar,
                        intcovar, weights, incl.markers,
                        maxit, tol, verbose, n.perm, n.node))
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
    temp <- scanone(cross, pheno.col=pheno.col, model=model, method=method, 
                    addcovar=addcovar, intcovar=intcovar, weights=weights,
                    maxit=maxit, tol=tol, verbose=FALSE)
    out.scanone <- temp[,-(1:2),drop=FALSE]
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
      if(chrtype[i]=="X" && (type=="bc" || type=="f2")) {
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
            reviseXdata(type, "full", sexpgm, geno=cross$geno[[i]]$data)
      }

    } # end loop over chromosomes
  } # end of if/else for method="mr" vs other 

  # columns in result matrix for each chromosome
  wh.col <- vector("list",n.chr)
  first.pos <- cumsum(c(1,n.pos))
  for(i in 1:n.chr)
    wh.col[[i]] <- seq(first.pos[i],by=1,length=n.pos[i])

  # initialize the results matrix
  if(n.phe > 1)
    results <- array(0,dim=c(sum(n.pos),sum(n.pos), n.phe))
  else
    results <- matrix(0,sum(n.pos),sum(n.pos))  

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
      if(i==j && chrtype[i]=="X") {
        col2drop <- dropXcol(type, sexpgm)
        n.col2drop <- sum(col2drop)
      }
      else {
        col2drop <- rep(0,n.gen[i]*n.gen[j])
        n.col2drop <- 0
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
                as.integer(n.phe),
                as.double(weights),
                result=as.double(rep(0,2*n.pos[i]*n.pos[j]*n.phe)),
                as.integer(n.col2drop),
                as.integer(col2drop),
                PACKAGE="qtl")
        z <- array(z$result,dim=c(n.pos[i], n.pos[j], 2*n.phe)) # rearrange the result 

        # update the final result matrix
        if(i == j) { # on same chromosome
          if(n.phe > 1)
            results[wh.col[[i]],wh.col[[j]],] <- z[,,1:n.phe]
          else
            results[wh.col[[i]],wh.col[[j]]] <- z[,,1]
          
          # do this just once: do null model and get neg log10 likelihood
          if(i==1) { 
            if(n.ac > 0)
              resid0 <- lm(pheno ~ ac, weights=weights^2)$resid
            else
              resid0 <- lm(pheno ~ 1, weights=weights^2)$resid
            sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
            nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
          }

#          cat("nllik0: ", nllik0, "\n")
        }
        else { # on different chromosomes
          if(n.phe > 1) {
            # full lod
            results[wh.col[[i]],wh.col[[j]],] <- z[,,1:n.phe]
            # epistasis lod - need to reshape the matrix
            results[wh.col[[j]],wh.col[[i]],] <- array(z[,,(n.phe+1):(2*n.phe)],
                                                       c(n.pos[j],n.pos[i], n.phe))
          }
          else { # only one phenotype, result is a matrix
            # full lod
            results[wh.col[[i]],wh.col[[j]]] <- z[,,1]
            # epistasis - need to reshape the matrix
            results[wh.col[[j]],wh.col[[i]]] <- matrix(z[,,2],n.pos[j],n.pos[i])
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
            if(method=="hk") {
              if(n.phe == 1)
                nllik0 <- (n.ind/2)*log10(sum((resid0*weights)^2))
              else # multiple phenotypes
                nllik0 <- apply(resid0, 2, function(x)
                                (n.ind/2)*log10(sum((x*weights)^2)))
            }
            else {
              sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
              nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
            }
          }


          if(verbose>1) cat("  --Calculating joint probs.\n")

          if(chrtype[i]=="X" && (type=="bc" || type=="f2")) {
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
          if(chrtype[i]=="X" && (type=="bc" || type=="f2")) 
            temp <- reviseXdata(type, "full", sexpgm, pairprob=temp)

          if(verbose>1) cat("  --Done.\n")

          if(method=="hk") {
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
                    as.integer(n.phe),
                    as.double(weights),
                    result=as.double(rep(0,n.pos[i]^2*n.phe)),
                    as.integer(n.col2drop),
                    as.integer(col2drop),
                    PACKAGE="qtl")
            ## fill results matrix
            if(n.phe == 1) 
              results[wh.col[[i]],wh.col[[i]]] <-
                matrix(z$result,ncol=n.pos[i])
            else # multiple phenotypes
              results[wh.col[[i]],wh.col[[i]],] <-
                array(z$result,c(n.pos[i],n.pos[i], n.phe))
          }
          else {
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
                    as.integer(n.col2drop),
                    as.integer(col2drop),
                    PACKAGE="qtl")
          # re-organize results
            results[wh.col[[i]],wh.col[[i]]] <-
              matrix(z$result,ncol=n.pos[i])
          }
          rm(temp) # remove the joint genotype probabilities
        } # end same chromosome
        
        else {
          if(method=="hk") {
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
                    as.integer(n.phe),
                    as.double(weights),
                    full=as.double(rep(0,n.pos[i]*n.pos[j]*n.phe)),
                    int=as.double(rep(0,n.pos[i]*n.pos[j]*n.phe)),
                    PACKAGE="qtl")
            ## reorgnize results
            if(n.phe == 1) {
              results[wh.col[[j]],wh.col[[i]]] <-
                matrix(z$full,ncol=n.pos[j])
              results[wh.col[[i]],wh.col[[j]]] <-
                matrix(z$int,ncol=n.pos[j])
            }
            else { # multiple phenotypes
              results[wh.col[[j]],wh.col[[i]],] <-
                array(z$full,c(n.pos[j], n.pos[i], n.phe))
              results[wh.col[[i]],wh.col[[j]],] <-
                array(z$int,c(n.pos[j], n.pos[i], n.phe))
            }
          }
          else  {
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
          }
        } # end different chromosome
      }
      else if(model=="binary" && method=="em") {
        if(i==j) { # same chromosome

          if(i==1) { # first time! do null model and get neg log10 likelihood
            if(n.ac > 0)
              nullfit <- glm(pheno ~ ac, family=binomial(link="logit"))
            else
              nullfit <- glm(pheno ~ 1, family=binomial(link="logit"))
            fitted <- nullfit$fitted
            nullcoef <- nullfit$coef
            nllik0 <- -sum(pheno*log10(fitted) + (1-pheno)*log10(1-fitted))
            if(verbose > 1) cat("null log lik: ", nllik0, "\n")
          }

          start <- c(rep(nullcoef[1],n.gen[i]),rep(0,n.gen[i]-1),
                     nullcoef[-1],rep(0,n.gen[i]*n.gen[i]+
                                      (n.gen[i]*(n.gen[i]-1)*n.ic)))

          if(n.col2drop)
            start <- c(start[!col2drop], rep(0,sum(col2drop)))

          if(verbose>1) cat("  --Calculating joint probs.\n")

          if(chrtype[i]=="X" && (type=="bc" || type=="f2")) {
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
          if(chrtype[i]=="X" && (type=="bc" || type=="f2")) 
            temp <- reviseXdata(type, "full", sexpgm, pairprob=temp)

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
                  as.integer(n.col2drop),
                  as.integer(col2drop),
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

        if(type=="f2") datai[datai>3] <- 0
        else if(type=="4way") datai[datai>4] <- 0

        if(chrtype[i]=="X" && (type=="bc" || type=="f2"))
            datai <- reviseXdata(type, "full", sexpgm, geno=datai)

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
                  as.integer(n.col2drop),
                  as.integer(col2drop),
                  PACKAGE="qtl")

          # re-organize results
          results[wh.col[[i]],wh.col[[i]]] <-
            matrix(z$result,ncol=n.pos[i])
        } # end same chromosome
        else {
          
          # replace missing and partially informative genotypes with 0's
          dataj <- cross$geno[[j]]$data
          dataj[is.na(dataj)] <- 0
          if(type=="f2") dataj[dataj>3] <- 0
          else if(type=="4way") dataj[dataj>4] <- 0

          if(chrtype[j]=="X" && (type=="bc" || type=="f2"))
            dataj <- reviseXdata(type, "full", sexpgm, geno=dataj)

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

  # subtract null neg log lik from lower tri
  if(method=="em")
    results[lower.tri(results)] <- nllik0 - results[lower.tri(results)]
  else if(method=="hk") {
    if(n.phe == 1)
      results[lower.tri(results)] <- nllik0 - results[lower.tri(results)]
    else { # multiple phenotypes
      lower <- lower.tri(results[,,1])
      for(itmp in 1:n.phe) {
        # I'm doing a loop here. I should put null model back to C function
        results[,,itmp][lower] <- nllik0[itmp] -
          results[,,itmp][lower]
      }
    }
  }

  # If the X chromosome was included, need to do an adjustment...
  scanoneX <- NULL
  if(any(gmap[,4])) { # the X chromosome was included

    # determine which covariates belong in null hypothesis
    temp <- scanoneXnull(type, sexpgm)
    adjustX <- temp$adjustX
    dfX <- temp$dfX
    sexpgmcovar <- temp$sexpgmcovar
      
    if(adjustX) {
      if(method=="mr" && any(is.na(pull.geno(cross)))) 
        warning("Scantwo with the X chr doesn't work quite right when method=\"mr\"\n",
                "    when there is missing genotype data.")

      if(model=="binary") {
        if(n.ac > 0)
          nullfitX <- glm(pheno ~ ac+sexpgmcovar, family=binomial(link="logit"))
        else
          nullfitX <- glm(pheno ~ sexpgmcovar, family=binomial(link="logit"))
        fittedX <- nullfitX$fitted
        nullcoefX <- nullfitX$coef
        nllikX <- -sum(pheno*log10(fittedX) + (1-pheno)*log10(1-fittedX))
        if(verbose > 1) cat("X chr null log lik: ", nllikX, "\n")
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

        if(method=="hk") {
          if(n.phe==1) 
            nllikX <- (n.ind/2)*log10(sum((residX*weights)^2))
          else 
            nllikX <- (n.ind/2)*apply(residX, 2, function(a,b)
                                    log10(sum((a*b)^2)), weights)
        }
        else {
          if(method=="imp" || method=="mr") {
            if(n.ac > 0) {
              out0 <- lm(pheno ~ ac, weights=weights^2)
              resid0 <- out0$resid
            }
            else {
              out0 <- lm(pheno ~ 1, weights=weights^2)
              resid0 <- out0$resid
            }
          
            if(n.phe > 1) {
              sig0 <- sqrt(apply(resid0, 2, function(a,b) sum((a*b)^2),weights)/n.ind)
              nllik0 <- sig0
              for(i in seq(along=nllik0))
                nllik0[i] <- -sum(dnorm(resid0[,i],0,sig0[i]/weights,log=TRUE))/log(10)
            }
            else {
              sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
              nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
            }
          }
 
          if(n.phe > 1) {
            sigX <- sqrt(apply(residX, 2, function(a,b) sum((a*b)^2),weights)/n.ind)
            nllikX <- sigX
            for(i in seq(along=nllikX))
              nllikX[i] <- -sum(dnorm(residX[,i],0,sigX[i]/weights,log=TRUE))/log(10)
          }
          else {
            sigX <- sqrt(sum((residX*weights)^2)/n.ind)
            nllikX <- -sum(dnorm(residX,0,sigX/weights,log=TRUE))/log(10)
          }
        }
      }

      if(n.phe > 1)
        wh <- ((gmap[row(results[,,1]),4] | gmap[col(results[,,1]),4]) &
               lower.tri(results[,,1]))
      else 
        wh <- ((gmap[row(results),4] | gmap[col(results),4]) & lower.tri(results))

      if(n.phe > 1) {
        for(i in 1:n.phe) 
          results[,,i][wh] <- results[,,i][wh] + nllikX[i] - nllik0[i]
      }
      else 
        results[wh] <- results[wh] + nllikX - nllik0

      if(run.scanone) {
        notxchr <- which(sapply(cross$geno,class)!="X")
        if(length(notxchr) > 0) {
          if(verbose) cat(" --Running scanone with special X chr covariates\n")
          temp <- scanone(subset(cross,chr=notxchr),
                          pheno.col=pheno.col, 
                          model=model, method=method,
                          addcovar=cbind(ac,sexpgmcovar),
                          intcovar=ic, weights=weights,
                          maxit=maxit, tol=tol, verbose=FALSE)

          scanoneX <- temp[,-(1:2),drop=FALSE]

          scanoneX <- rbind(scanoneX,
                            out.scanone[rownames(gmap),,drop=FALSE][gmap[,4],,drop=FALSE])

          scanoneX <- scanoneX[rownames(gmap),,drop=FALSE]
        }
        else {
          scanoneX <- out.scanone[rownames(gmap),,drop=FALSE][gmap[,4],,drop=FALSE]
          scanoneX <- scanoneX[rownames(gmap),,drop=FALSE]
        }
      }
    }
  }
  
#  if(any(is.na(results) | results < -1e-6 | results == Inf))
#    warning("Some LOD scores NA, Inf or < 0")
  if(any(is.na(results)))
    warning("Some LOD scores NA")
  if(any(!is.na(results) & results < 0))
    warning("Some LOD scores < 0")
  if(any(!is.na(results) & (results == Inf | results == -Inf)))
    warning("Some LOD scores = Inf or -Inf")
  
  if(!is.null(scanoneX)) scanoneX <- as.matrix(scanoneX)

  # output has 2 fields, lod and map
  out <- list(lod=results,map=gmap,scanoneX=scanoneX)
  class(out) <- "scantwo"

  # fill in scanone result
  if(run.scanone) { 
    if(n.phe == 1)
      diag(out$lod) <- out.scanone[rownames(out$map),]
    else {
      for(iphe in 1:n.phe) 
        diag(out$lod[,,iphe]) <- out.scanone[rownames(out$map),iphe]
    }
  }

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
function(cross, pheno.col=1, model=c("normal","binary"),
         method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
         addcovar=NULL, intcovar=NULL, weights=NULL,
         incl.markers=FALSE, maxit=4000, tol=1e-4, verbose=FALSE,
         n.perm=1000, n.node=1) 
{
  method <- match.arg(method)
  model <- match.arg(model)

  # for running on clusters, check if the setup is correct
  if(n.node > 1) { ## run on cluster
    ## initialize cluster
    if(!require(snow))
      stop("You have to install SNOW package to use cluster")
    ## make cluster object
    cl <- makeMPIcluster(nnodes)
    ## correct the possible correlation problem
    clusterApply(cl, runif(length(cl), max=1000000000), set.seed)
    ## turn off verbose
    verbose <- FALSE
  }
  
  ## start permutation
  if(n.node > 1) { ## run on clusters
    ## calculate the number of permutation needed in each node
    nperm.cluster <- rep(floor((n.perm-1)/n.node), n.node)
    ## maybe some leftovers
    leftover <- n.perm - 1 - sum(nperm.cluster)
    if(leftover > 0)
      nperm.cluster[1:leftover] <- nperm.cluster[1:leftover] + 1
    ## load library on all nodes
    clusterEvalQ(cl, library(qtl))
    ## use clusterApply to run permutation on all nodes
    ## note that the only different parameter passed to function
    ## is nperm.cluster.
    cat(paste("Doing permutation on", nnodes, "cluster nodes ... \n"))
    pstar.nodes <- clusterApply(cl, nperm.cluster, scantwo.perm.engine,
                                cross, pheno.col, model, method,
                                addcovar, intcovar, weights,
                                incl.markers, maxit, tol, verbose)
    
    ## collect the result from clusters
    ## this part is to be done !!! - hao wu
  }
  else { ## on single node
    res <- scantwo.perm.engine(n.perm, cross, pheno.col, model, 
                               method, addcovar, intcovar, weights,
                               incl.markers, maxit, tol, verbose)
  }
  
  res
}


######################################################################
#
# Engine function to scantwo permutation
#
######################################################################
scantwo.perm.engine <-
function(n.perm, cross, pheno.col, model,
         method, addcovar, intcovar, weights,
         incl.markers, maxit, tol, verbose)
{
  ## local variables
  n.phe <- length(pheno.col)
#  n.addcov <- ncol(addcovar)
  n.ind <- dim(cross$pheno)[1]

  ## initialize result
  ## it's a list of two components, LOD.jnt and LOD.interxn
  LOD.jnt <- matrix(0, n.perm, n.phe)
  LOD.interxn <- matrix(0, n.perm, n.phe)
#  res <- matrix(ncol=2,nrow=n.perm)
      
  if( (n.phe == 1) && ((method=="imp") || (method=="hk")) &&
       is.null(addcovar) && is.null(intcovar) ) {
    if(verbose)
      cat("Doing permutation in batch mode ...\n")
    ## if there's only one phenotype, no covariate, and method is imp or hk,
    ## generate permuted phenotype as a matrix and do permutation
    ## as multiple phenotypes
    ord <- matrix(0, n.ind, n.perm)
    for(iperm in 1:n.perm)
      ord[,iperm] <- sample(n.ind)

    cross$pheno <- matrix(cross$pheno[,pheno.col][ord], nrow=n.ind)

    pheno.col <- 1:n.perm
    tem <- scantwo(cross,,pheno.col,model,method,addcovar,
                   intcovar,weights,run.scanone=FALSE,
                   maxit=maxit,incl.markers=incl.markers,
                   tol=tol, verbose=TRUE, n.perm=-1)
    ## find the maximum LOD on each permutation
    for(iperm in 1:n.perm) {
      tmplod <- tem$lod[,,iperm]
      LOD.jnt[iperm, 1] <- max(tmplod[tmplod < Inf & row(tmplod)>col(tmplod)],
                               na.rm=TRUE )
      LOD.interxn[iperm,1] <- max(tmplod[tmplod < Inf & row(tmplod)<col(tmplod)],
                                  na.rm=TRUE )
    }
  }
  else { ## all other cases, do one permutation at a time
    if(method=="mr-imp") # save version with missing genotypes
      tempcross <- cross
    if(method=="mr-argmax") # impute genotypes
      cross <- fill.geno(cross,method="argmax")
    addcovarp <- intcovarp <- NULL
    if(!is.null(addcovar)) addcovar <- as.matrix(addcovar)
    if(!is.null(intcovar)) intcovar <- as.matrix(intcovar)
    ## permutation loop
    for(i in 1:n.perm) {
      if(verbose) cat("Permutation", i, "\n")
      ## impute genotypes for method "mr-imp"
      if(method=="mr-imp") cross <- fill.geno(tempcross)

      o <- sample(1:n.ind)
      cross$pheno <- cross$pheno[o,,drop=FALSE]
      if(!is.null(addcovar)) addcovarp <- addcovar[o,,drop=FALSE]
      if(!is.null(intcovar)) intcovarp <- intcovar[o,,drop=FALSE]
      tem <- scantwo(cross,  pheno.col=pheno.col, model=model, 
                     method=method, addcovar=addcovarp,
                     intcovar=intcovarp, incl.markers=incl.markers,
                     weights=weights, run.scanone=FALSE, maxit=maxit,
                     tol=tol, verbose=FALSE, n.perm = -1)

      ## calculate max joint lod and interactive lod for each phenotype
      ## it's to take max of the two triangles for each trait
      for(iphe in 1:n.phe) {
        if(n.phe == 1)
          tmplod <- tem$lod
        else
          tmplod <- tem$lod[,,iphe]
        LOD.jnt[i, iphe] <- max(tmplod[tmplod < Inf & row(tmplod)>col(tmplod)],
                                na.rm=TRUE )
        LOD.interxn[i,iphe] <- max(tmplod[tmplod < Inf & row(tmplod)<col(tmplod)],
                                   na.rm=TRUE )
      }
    }

  }

  ## make result
  res <- list(LOD.jnt=LOD.jnt, LOD.interxn=LOD.interxn)
  attr(res,"method") <- method
  # return
  res
}

######################################################################
#
# summerize the result from scantwo
#
######################################################################
summary.scantwo <-
function (object, thresholds = c(0, 0, 0), lodcolumn=1,
          type = c("joint","interaction"), ...)
{
  type <- match.arg(type)
  
  if(length(dim(object$lod)) > 2) { # results from multiple phenotypes
    if(length(lodcolumn) > 1) {
      warning("Argument lodcolumn should be of length 1.")
      lodcolumn <- lodcolumn[1]
    }
      
    if(lodcolumn < 0 || lodcolumn > dim(object$lod)[3])
      stop("Argument lodcolumn misspecified.")
    object$lod <- object$lod[,,lodcolumn]
  }

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
      tmplod1 <- lod[wh.index[[j]], wh.index[[i]],drop=FALSE]
      if(!is.null(condlod)) {
        if(i==j) tmpcondlod <- condlod[wh.index[[i]],wh.index[[i]],drop=FALSE]
        else {
          tmpcondlod1 <- condlod[wh.index[[j]],wh.index[[i]],drop=FALSE]
          tmpcondlod2 <- condlod[wh.index[[i]],wh.index[[j]],drop=FALSE]
        }
      }

      if(i != j) tmplod2 <- lod[wh.index[[i]], wh.index[[j]],drop=FALSE]
      else tmplod2 <- tmplod1
      

      if(type == "joint") {
        if(i == j) {
          tri <- lower.tri(tmplod1)
          lod.joint <- max(tmplod1[tri])
          idx <- which(tmplod1 == lod.joint & tri, arr.ind=TRUE)
          if(!is.matrix(idx)) {
            cat("problem\n")
            return(tmplod1)
          }
        }
        else {
          lod.joint <- max(tmplod1)
          idx <- which(tmplod1 == lod.joint, arr.ind=TRUE)
          if(!is.matrix(idx)) {
            cat("problem\n")
            return(tmplod1)
          }
        }
        if(nrow(idx)>1) idx <- idx[sample(nrow(idx),1),,drop=FALSE]
        idx.row <- idx[1]
        idx.col <- idx[2]
        
        lod.int <- tmplod2[idx.col, idx.row,drop=FALSE]
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
        if(nrow(idx)>1) idx <- idx[sample(nrow(idx),1),,drop=FALSE]
        idx.row <- idx[2]
        idx.col <- idx[1]
        
        lod.joint <- tmplod1[idx.row, idx.col,drop=FALSE]
      }
      
      full.idx.row <- idx.row + wh.index[[j]][1] - 1
      full.idx.col <- idx.col + wh.index[[i]][1] - 1

      flag <- FALSE # a flag to indicate whether there's any peak on this pair
      if(lod.joint >= thrfull) {
        if(includes.scanone) {
          if(i==j) {
            lod.q1 <- tmpcondlod[idx.row,idx.col,drop=FALSE]
            lod.q2 <- tmpcondlod[idx.col,idx.row,drop=FALSE]
          }
          else {
            lod.q1 <- tmpcondlod1[idx.row,idx.col,drop=FALSE]
            lod.q2 <- tmpcondlod2[idx.col,idx.row,drop=FALSE]
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
