######################################################################
#
# scantwo.R
#
# copyright (c) 2001-2, Karl W Broman, Johns Hopkins University,
#                       and Hao Wu, The Jackson Lab.
# last modified June, 2002
# first written Nov, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Hao Wu (The Jackson Lab) wrote the initial code for the imputation
# method and the plot.scantwo and summary.scantwo functions.
#
# Part of the R/qtl package
# Contains: scantwo, plot.scantwo, scantwo.perm, summary.scantwo
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
         addcovar=NULL, intcovar=NULL, x.treatment=c("simple","full"),
         run.scanone=TRUE, incl.markers=FALSE, maxit=4000, tol=1e-4,
         trace=TRUE, n.perm)
{
  if(method=="im") # warning in case old terminology is used
    warning("Method \"im\" is now called \"em\"; running method \"imp\".")
  method <- match.arg(method)
  
  origcross <- cross

  # pull out chromosomes to be scanned
  if(!missing(chr)) cross <- subset(cross,chr=chr)

  # check phenotypes and covariates; drop individuals with missing values
  # in case of permutation test, only do checks once
  if(missing(n.perm) || n.perm>0) { 
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

  # Problems with EX w/ X chromosome: just use H-K for now.
  if(any(chrtype=="X") && method=="em") {
    sexpgm <- getsex(cross)
    if(!is.null(sexpgm$sex) || !is.null(sexpgm$pgm)) {
      warning("EM not working for X chromosomes; using H-K instead.")
      method <- "hk"
    }
  }

  if(missing(n.perm) || n.perm == 0) { # not in the midst of permutations
    if(method=="mr-argmax")
      cross <- fill.geno(cross,method="argmax")
    if(method=="mr-imp")
      cross <- fill.geno(cross,method="imp")
  }

  # if n.perm specified, do a permutation test
  if(!missing(n.perm) && n.perm>0) { 
    return(scantwo.perm(cross, pheno.col, method, addcovar,
                        intcovar, x.treatment, incl.markers, maxit, tol,
                        trace, n.perm))
  }

  if(run.scanone) { # also do scanone
    if(trace) cat(" --Running scanone\n")
    temp <- scanone(cross,, pheno.col, , method, addcovar, intcovar,
                    x.treatment, maxit=maxit, tol=tol, trace=FALSE)
    nam <- rownames(temp)
    out.scanone <- temp[,3]
    names(out.scanone) <- nam
    if(trace) cat(" --Running scantwo\n")
  }

  if(method=="mr" || method=="mr-imp" || method=="mr-argmax") { # marker regression
    # number of genotypes on each chromosome, 
    #     combine the genetic maps for all chromosomes
    map <- unlist(pull.map(cross))
    names(map) <- unlist(lapply(pull.map(cross),names))
    n.pos <- nmar(cross)
    gmap <- data.frame(chr=rep(names(cross$geno),n.pos),
                       pos=map,
                       eq.spacing=rep(1,sum(n.pos)))

    # number of possible genotypes for each chromosome
    n.gen <- 1:n.chr
    for(i in 1:n.chr) { 
      if(chrtype[i]=="X") sexpgm <- getsex(cross)
      else sexpgm <- NULL

      gen.names <- getgenonames(type, chrtype[i], x.treatment, sexpgm)
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
      if(chrtype[i]=="X") sexpgm <- getsex(cross)
      else sexpgm <- NULL

      gen.names <- getgenonames(type, chrtype[i], x.treatment, sexpgm)
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
      if(length(o) > 0) 
        w[o] <- paste(w[o],names(cross$geno)[i],sep=".c")
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
      gmap <- rbind(gmap, cbind(map,eq.spacing=eq.sp.pos))
      n.pos[i] <- length(keep.pos[[i]])

      # Fix genotype probabilities or imputations for X chromosome
      if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss")) {
        if(method=="imp") 
          cross$geno[[i]]$draws <-
            fixXdata(type, x.treatment, sexpgm, draws=cross$geno[[i]]$draws)
        else
          cross$geno[[i]]$prob <-
            fixXdata(type, x.treatment, sexpgm, prob=cross$geno[[i]]$prob)
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

      # print the current working pair
      if(trace) cat(paste(" (", names(cross$geno)[i], ",",
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
                as.double(addcovar),
                as.integer(n.addcovar),
                as.double(intcovar),
                as.integer(n.intcovar),
                as.double(pheno),
                result=as.double(rep(0,2*n.pos[i]*n.pos[j])),
                PACKAGE="qtl")
        z <- array(z$result,dim=c(n.pos[i], n.pos[j], 2)) # rearrange the result 

        # update the final result matrix
        results[wh.col[[i]],wh.col[[j]]] <- z[,,1]
        if(i != j) results[wh.col[[j]],wh.col[[i]]] <- t(z[,,2])
      }
      else if(method=="hk" || method=="em") {
        if(i==j) { # same chromosome

          if(i==1) { # first time! do null model and get neg log10 likelihood
            if(n.addcovar > 0) resid0 <- lm(pheno ~ addcovar)$resid
            else resid0 <- pheno - mean(pheno)
            if(method=="hk") nllik0 <- (n.ind/2)*log10(sum(resid0^2))
            else {
              sig0 <- sqrt(sum(resid0^2)/n.ind)
              nllik0 <- -sum(dnorm(resid0,0,sig0,log=TRUE))/log(10)
            }
          }

          if(trace>1) cat("  --Calculating joint probs.\n")
          # calculate joint genotype probabilities for all pairs of positions
          stp <- attr(cross$geno[[i]]$prob, "step")
          oe <- attr(cross$geno[[i]]$prob, "off.end")
          err <- attr(cross$geno[[i]]$prob, "error.prob")
          mf <- attr(cross$geno[[i]]$prob, "map.function")
          temp <- calc.pairprob(subset.cross(cross,chr=i),stp,oe,err,mf)

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

          if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss")) {
            temp <- fixXdata(type, x.treatment, sexpgm, pairprob=temp)
            temp[temp==0] <- 1e-5 # << temp fix for problems with X chromosome
          }

          if(trace>1) cat("  --Done.\n")

          if(method=="hk") 
            z <- .C("R_scantwo_1chr_hk", 
                    as.integer(n.ind),
                    as.integer(n.pos[i]),
                    as.integer(n.gen[i]),
                    as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                    as.double(temp),
                    as.double(addcovar),
                    as.integer(n.addcovar),
                    as.double(intcovar),
                    as.integer(n.intcovar),
                    as.double(pheno),
                    result=as.double(rep(0,n.pos[i]^2)),
                    PACKAGE="qtl")
          else
            z <- .C("R_scantwo_1chr_em", 
                    as.integer(n.ind),
                    as.integer(n.pos[i]),
                    as.integer(n.gen[i]),
                    as.double(temp),
                    as.double(addcovar),
                    as.integer(n.addcovar),
                    as.double(intcovar),
                    as.integer(n.intcovar),
                    as.double(pheno),
                    result=as.double(rep(0,n.pos[i]^2)),
                    as.integer(maxit),
                    as.double(tol),
                    as.integer(trace),
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
                    as.double(addcovar),
                    as.integer(n.addcovar),
                    as.double(intcovar),
                    as.integer(n.intcovar),
                    as.double(pheno),
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
                    as.double(addcovar),
                    as.integer(n.addcovar),
                    as.double(intcovar),
                    as.integer(n.intcovar),
                    as.double(pheno),
                    full=as.double(rep(0,n.pos[i]*n.pos[j])),
                    int=as.double(rep(0,n.pos[i]*n.pos[j])),
                    as.integer(maxit),
                    as.double(tol),
                    as.integer(trace),
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

        if(chrtype[i]=="X" && (type=="bc" || type=="f2" || type=="f2ss"))
          datai <- fixXdata(type, x.treatment, sexpgm, geno=datai)

        if(i==j) { # same chromosome

          z <- .C("R_scantwo_1chr_mr",
                  as.integer(n.ind),
                  as.integer(n.pos[i]),
                  as.integer(n.gen[i]),
                  as.integer(datai),
                  as.double(addcovar),
                  as.integer(n.addcovar),
                  as.double(intcovar),
                  as.integer(n.intcovar),
                  as.double(pheno),
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

          if(chrtype[j]=="X" && (type=="bc" || type=="f2" || type=="f2ss"))
            dataj <- fixXdata(type, x.treatment, sexpgm, geno=dataj)

          z <- .C("R_scantwo_2chr_mr",
                  as.integer(n.ind),
                  as.integer(n.pos[i]),
                  as.integer(n.pos[j]),
                  as.integer(n.gen[i]),
                  as.integer(n.gen[j]),
                  as.integer(datai),
                  as.integer(dataj),
                  as.double(addcovar),
                  as.integer(n.addcovar),
                  as.double(intcovar),
                  as.integer(n.intcovar),
                  as.double(pheno),
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

  if(any(is.na(results) | results < -1e-6 | results == Inf))
    warning("Some LOD scores NA, Inf or < 0")
#  results[is.na(results) | results<0 | results == Inf] <- 0
  

  # output has 2 fields, lod and map
  out <- list(lod=results,map=gmap)
  class(out) <- "scantwo"

  if(run.scanone) { # also do scanone
    if(method=="mr" || method=="mr-imp" || method=="mr-argmax")
      diag(out$lod) <- out.scanone
    else diag(out$lod) <- out.scanone[rownames(out$map)]
  } # end scanone 

  attr(out,"method") <- method
  attr(out,"type") <- type
  out
}

######################################################################
#
# plot.scantwo: plot output from scantwo
#
######################################################################

plot.scantwo <- 
function(x,chr,incl.markers=FALSE,zlim,contours=FALSE,
         main,zscale=TRUE,...)
{
#  if( !any(class(x) == "scantwo") )
#    stop("Input variable is not an object of class scantwo!")
  
  lod <- x$lod
  map <- x$map

  # deal with bad LOD score values
  if(any(is.na(lod) | lod< -1e-6 | lod==Inf)) 
    warning("Some LOD scores NA, Inf or < 0; set to 0")
  lod[is.na(lod) | lod<0 | lod == Inf] <- 0

  # if incl.markers is FALSE, drop positions
  #     for which third column of map is 0
  if(!incl.markers && any(map[,3]==0)) {
    o <- (map[,3]==1)
    lod <- lod[o,o]
    map <- map[o,]
  }

  # turn NA's into negative values (plotted as white?)
  lod[is.na(lod)] <- -10

  # pull out desired chromosomes
  if(missing(chr) || length(chr)==0)
    chr <- unique(map[,1])
  else {
    a <- unique(map[,1])
    if(is.numeric(chr) && all(chr < 0)) 
      chr <- a[-match(-chr,a)]
    else chr <- a[match(chr,a)]
    keep <- (1:nrow(map))[!is.na(match(map[,1],chr))]

    map <- map[keep,]
    lod <- lod[keep,keep]
  }

  chr <- as.character(chr)

  if( missing(zlim) ) { # no given zlim
    # calculate the zlim for interactive and joint
    zlim.int <- max( lod[row(lod)<col(lod)] )
    zlim.jnt <- max( lod[row(lod)>=col(lod)] )
  }
  else {
    zlim.int <- zlim[2]
    zlim.jnt <- zlim[1]
  }

  
  # rescale the data in upper triangle based on zlims.jnt
  lod[row(lod)<col(lod)] <- lod[row(lod)<col(lod)]*zlim.jnt/zlim.int

  if(missing(zlim)) zlim.jnt <- max(lod)

  # make sure LOD values are below (0,zlim.jnt) or update zlim.jnt
  if(max(lod) > zlim.jnt) {
    warning("LOD values out of range; updating zlim.")
    temp <- max(lod)
    zlim.int <- zlim.int*temp/zlim.jnt
    zlim.jnt <- temp
  }

  # save old par parameters, to restore them on exit
  old.mar <- par("mar")
  old.las <- par("las")
  old.mfrow <- par("mfrow") 
  on.exit(par(las=old.las,mar=old.mar,mfrow=old.mfrow))
  par(las=1)

  if(zscale) {
    layout(cbind(1,2),c(6,1))
    par(mar=c(5,4,4,2)+0.1)
  }

  image( 1:ncol(lod), 1:nrow(lod), lod, ylab="Positions", xlab="Positions",
         zlim=c(0,zlim.jnt), col=rev(rainbow(256,start=0,end=2/3)) )

  # add contours if requested
  if(contours) contour(1:ncol(lod), 1:nrow(lod), lod, add=TRUE)

  # calculate how many markers in each chromesome
  n.mar <- NULL
  for ( i in 1:length(chr) )
    n.mar[i] <- sum(map[,1]==chr[i])
  
  # plot lines at the chromosome boundaries
  wh <- c(0.5,cumsum(n.mar)+0.5)
  abline(v=wh,xpd=FALSE)
  abline(h=wh,xpd=FALSE)

  # add chromesome numbers
  a <- par("usr")
  for(i in 1:length(n.mar)) 
    text(mean(wh[i+c(0,1)]),a[4]+(a[4]-a[3])*0.025,chr[i],xpd=TRUE)
  for(i in 1:length(n.mar)) 
    text(a[2]+(a[2]-a[1])*0.025,mean(wh[i+c(0,1)]),chr[i],xpd=TRUE)
  
  # add title
  if(!missing(main)) title(main=main)

  if(zscale) {
    # plot the colormap
    par(mar=c(5,2,4,2)+0.1)
    colorstep <- zlim.jnt/255
    image( x=1:1, y=seq(0,zlim.jnt,colorstep), z=matrix(c(1:256),1,256),
          zlim=c(1,256),ylab="",xlab="", xaxt="n", yaxt="n",
          col=rev(rainbow(256,start=0,end=2/3)) )
    # make sure there's a box around it
    u <- par("usr") 
    abline(v=u[1:2],xpd=FALSE)
    abline(h=u[3:4],xpd=FALSE)
  
    # figure out how big the axis labels should be
    fin <- par("fin")[1] # figure width in inches
    pin <- par("pin")[1] # plot width in inches
    mai <- par("mai")[2] # margin width in inches
                         # note: pin + 2*mai = fin
    xlen.mar <- mai/pin*diff(u[1:2])

    # axis for joint LODs
    yloc <- pretty(c(0,zlim.jnt),4)
    yloc <- yloc[yloc<=u[4]]
    segments(u[2],yloc,u[2]+xlen.mar/4,yloc,xpd=TRUE)
    text(u[2]+xlen.mar/3,yloc,as.character(yloc),xpd=TRUE,adj=0)
 
    # axis for int've LODs
    yloc <- pretty(c(0,zlim.int),4)
    yloc.rev <- yloc*zlim.jnt/zlim.int
    yloc <- yloc[yloc.rev <= u[4]]
    yloc.rev <- yloc.rev[yloc.rev<=u[4]]
    segments(u[1],yloc.rev,u[1]-xlen.mar/4,yloc.rev,xpd=TRUE)
    text(u[1]-xlen.mar/3,yloc.rev,as.character(yloc),xpd=TRUE,adj=1)
  }
}

######################################################################
#
# scantwo.perm: Permutation test of scantwo
#
######################################################################

scantwo.perm <-
function(cross, pheno.col=1,
         method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
         addcovar=NULL, intcovar=NULL, x.treatment=c("simple","full"),
         incl.markers=FALSE, maxit=4000, tol=1e-4, trace=FALSE,
         n.perm=1000) 
{
  method <- match.arg(method)
  x.treatment <- match.arg(x.treatment)

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
    if(trace) cat("Permutation", i, "\n")

    # impute genotypes for method "mr-imp"
    if(method=="mr-imp") cross <- fill.geno(tempcross)

    o <- sample(1:n.ind)
    cross$pheno <- cross$pheno[o,,drop=FALSE]
    if(!is.null(addcovar)) addcovarp <- addcovar[o,,drop=FALSE]
    if(!is.null(intcovar)) intcovarp <- intcovar[o,,drop=FALSE]
    tem <- scantwo(cross,  pheno.col=pheno.col,
                   method=method, addcovar=addcovarp,
                   intcovar=intcovarp, incl.markers=incl.markers,
                   x.treatment=x.treatment,
                   run.scanone=FALSE, maxit=maxit, tol=tol,
                   trace=FALSE, n.perm = -1)

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
function(object, thresholds=c(0,0,0), ...)
{
  if(length(thresholds) < 3) {
    if(length(thresholds)==1) thresholds <- c(thresholds,0,0)
    else stop("You must give three thresholds: full, interaction and main\n")
  }
    
  thrfull <- thresholds[1]
  thrint <- thresholds[2]
  thrmain <- thresholds[3]

  lod <- object$lod
  map <- object$map

  # deal with bad LOD score values
  if(any(is.na(lod) | lod < -1e-6 | lod==Inf)) 
    warning("Some LOD scores NA, Inf or < 0; set to 0")
  lod[is.na(lod) | lod<0 | lod == Inf] <- 0

  # if there's no mainscan result, ignore the thresholds
  #     and don't include the 4 conditional LOD columns
  if(all(is.na(diag(lod)) | diag(lod) < 1e-10)) includes.scanone <- FALSE
  else includes.scanone <- TRUE

  crosstype <- attr(object,"type")
  if(is.null(crosstype)) {
    warning("No type attribute in input data; assuming backcross.")
    crosstype <- "bc"
  }

  # calculate the degree of freedom
  if(crosstype == "bc" || crosstype=="riself" || crosstype=="risib") {
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
  else stop(paste("Don't know what to do with cross type", crosstype))
  
  # chromsomes in the result
  chr <- unique(map[,1])
  n.chr <- length(chr)
  
  # calculate the locations of each chromosome within the LOD matrix
  wh.index <- vector("list",n.chr)
  n <- nrow(map)
  for(i in 1:n.chr)
    wh.index[[i]] <- (1:n)[map[,1]==chr[i]]
    
  results <- NULL
  # go through each pair of chromosomes
  for (i in 1:n.chr) {
    for (j in i:n.chr) {

      if( i == j) {# same chromesome
        tmplod <- lod[wh.index[[i]],wh.index[[i]]]
        # find the maximum LOD score and its location [idx.row, idx.col]
        lod.joint <- max( tmplod )
        idx.row <- row(tmplod)[tmplod==lod.joint][1]
        idx.col <- col(tmplod)[tmplod==lod.joint][1]
        # get the LOD score for H{f} vs. H{1+2}
        lod.int <- tmplod[ idx.col,idx.row ]
      }
      else { # different chromesome
        tmplod1 <- lod[wh.index[[j]], wh.index[[i]]]
        tmplod2 <- lod[wh.index[[i]], wh.index[[j]]]
        # find the maximum LOD score and its location [idx.row, idx.col]
        lod.joint <- max(tmplod1)
        idx.row <- row(tmplod1)[tmplod1==lod.joint][1]
        idx.col <- col(tmplod1)[tmplod1==lod.joint][1]
        # get the LOD score for H{f} vs. H{1+2}
        lod.int <- tmplod2[ idx.col,idx.row ]
      }

      # convert back to position in big matrix
      idx.row <- idx.row + wh.index[[j]][1] - 1
      idx.col <- idx.col + wh.index[[i]][1] - 1
      
      if(lod.joint >= thrfull) {
        if(includes.scanone) {
          # lod(Q1|Q2) comparing Q1+Q2 (no intxn) to Q2 alone
          lod.q1 <- lod.joint - lod.int - lod[idx.row,idx.row]
          # lod(Q2|Q1) comparing Q1+Q2 (no intxn) to Q1 alone
          lod.q2 <- lod.joint - lod.int - lod[idx.col,idx.col]
        
          if(lod.int >= thrint || min(lod.q1,lod.q2) >= thrmain) {
            i.pos <- map[idx.col,2]
            j.pos <- map[idx.row,2]

            results <-
              rbind(results,
                    data.frame(chr[i], chr[j], i.pos, j.pos, lod.joint,
                               1-pchisq(2*log(10)*lod.joint, df.int+2*df.add),
                               lod.int, 1-pchisq(2*log(10)*lod.int,df.int),
                               lod.q1, 1-pchisq(2*log(10)*lod.q1,df.add),
                               lod.q2, 1-pchisq(2*log(10)*lod.q2,df.add)))
          }
        }
        else { # no scanone available
          i.pos <- map[idx.col,2]
          j.pos <- map[idx.row,2]
          results <- rbind(results,
                           data.frame(chr[i], chr[j], i.pos, j.pos, lod.joint,
                                      1-pchisq(2*log(10)*lod.joint, df.int+2*df.add),
                                      lod.int, 1-pchisq(2*log(10)*lod.int,df.int)))
        }
      }
    }
  }

  if(is.null(results)) {
    cat("    There were no pairs of loci meeting the criteria.\n")
    invisible()
  }
  else {
    if(includes.scanone)
      colnames(results) <- c("chr1","chr2","pos1","pos2",
                             "lod.joint","p.joint",
                             "lod.int","p.int","lod.q1","p.q1",
                             "lod.q2","p.q2")
    else
       colnames(results) <- c("chr1","chr2","pos1","pos2",
                              "lod.joint","p.joint",
                              "lod.int","p.int")

    results <- as.data.frame(results)
    class(results) <- c("summary.scantwo","data.frame")
    return(results)
  }
}

print.summary.scantwo <-
function(x,...)
{
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
