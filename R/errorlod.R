######################################################################
#
# errorlod.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# Oct, 2001; Sept, 2001; May, 2001; Apr, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: calc.errorlod, plot.errorlod, top.errorlod
#
######################################################################

######################################################################
#
# calc.errorlod: Calculate LOD scores indicating likely genotyping
#                errors.
#
######################################################################

calc.errorlod <-
function(cross, error.prob=0.01,
         map.function=c("haldane","kosambi","c-f"))
{

  # map function
  map.function <- match.arg(map.function)

  n.ind <- nind(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)

  if(class(cross)[1]=="bc") cfunc <- "calc_errorlod_bc"
  else if(class(cross)[1]=="f2") cfunc <- "calc_errorlod_f2"
  else if(class(cross)[1]=="4way") cfunc <- "calc_errorlod_4way"
  else stop(paste("calc.errorlod not available for cross type",
                  class(cross)[1],"."))
  
  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {

    # skip chromosomes with only 1 marker
    if(n.mar[i] < 2) next

    if(is.na(match("prob",names(cross$geno[[i]])))) {
      # need to run calc.genoprob
      warning("First running calc.genoprob.")
      cross <- calc.genoprob(cross,error.prob=error.prob,
                             map.function=map.function)
      Pr <- cross$geno[[i]]$prob
    }
    else {
      # if error.prob doesn't correspond to what was used when
      #     running calc.genoprob(), re-run calc.genoprob()
      if(abs(attr(cross$geno[[i]]$prob,"error.prob")
             - error.prob) > 1e-9) {
        warning("Re-running calc.genoprob()")
        cross <-
          calc.genoprob(cross,error.prob=error.prob,
                        step=attr(cross$geno[[i]]$prob,"step"),
                        off.end=attr(cross$geno[[i]]$prob,"off.end"),
                        map.function=map.function)
      }
         
      Pr <- cross$geno[[i]]$prob
      u <- grep("^loc\-*[0-9]+",colnames(Pr))
      if(length(u) > 0) Pr <- Pr[,-u,]
    }
    
    nm <- dim(Pr)[2]
    dat <- cross$geno[[i]]$data
    dat[is.na(dat)] <- 0
    
    z <- .C(cfunc,
            as.integer(n.ind),
            as.integer(nm),
            as.integer(dat),
            as.double(error.prob),
            as.double(Pr),
            errlod=as.double(rep(0,n.ind*nm)),
            PACKAGE="qtl")

    errlod <- array(z$errlod, dim=dim(Pr)[1:2])

    dimnames(errlod) <- list(NULL,colnames(cross$geno[[i]]$data))
    cross$geno[[i]]$errorlod <- errlod

    # attribute set to the error.prob value used, for later
    #     reference.
    attr(cross$geno[[i]]$errorlod,"error.prob") <- error.prob
  }

  cross
}

  

  

######################################################################
#
# plot.errorlod
#
######################################################################

plot.errorlod <-
function(x, chr, ...)
{
  cross <- x
  if(!missing(chr)) cross <- pull.chr(cross,chr)

  # remove chromosomes with < 2 markers
  n.mar <- nmar(cross)
  cross <- pull.chr(cross,names(n.mar)[n.mar >= 2])
  n.chr <- nchr(cross)

  errlod <- NULL
  for(i in 1:n.chr) {
    if(is.na(match("errorlod",names(cross$geno[[i]])))) { # need to run calc.errorlod
      warning("First running calc.errorlod.")
      cross <- calc.errorlod(cross,error.prob=0.01,map.function="haldane")
    }
    errlod <- cbind(errlod,cross$geno[[i]]$errorlod)
  }

  errlod <- t(errlod)

  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=TRUE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))

  colors <- c("white", "gray85", "hotpink", "firebrick")

  # plot grid 
  image(1:nrow(errlod),1:ncol(errlod),errlod,
        ylab="Individuals",xlab="Markers",col=colors,
        breaks=c(-1,1,2,3,max(errlod)))

  # plot lines at the chromosome boundaries
  n.mar <- nmar(cross)
  chr.names <- names(cross$geno)
  a <- c(0.5,cumsum(n.mar)+0.5)

  # the following makes the lines go slightly above the plotting region
  b <- par("usr")
  segments(a,b[3],a,b[4]+diff(b[3:4])*0.02)

  # this line adds a line above and below the image
  #     (the image function seems to leave these out)
  abline(h=0.5+c(0,ncol(errlod)),xpd=FALSE)

  # add chromosome numbers
  a <- par("usr")
  wh <- cumsum(c(0.5,n.mar))
  for(i in 1:length(n.mar)) 
    text(mean(wh[i+c(0,1)]),a[4]+(a[4]-a[3])*0.025,chr.names[i])

  title(main="Genotyping error LOD scores")

}


######################################################################
#
# top.errorlod
#
# Picks out the genotypes having errorlod values above some cutoff
#
######################################################################

top.errorlod <-
function(cross, chr, cutoff=2, msg=TRUE)  
{
  if(!missing(chr)) cross <- pull.chr(cross,chr)

  mar <- ind <- lod <- chr <- NULL

  # remove chromosomes with < 2 markers
  n.mar <- nmar(cross)
  cross <- pull.chr(cross,names(n.mar)[n.mar >= 2])

  flag <- 0
  for(i in 1:nchr(cross)) {
    
    if(is.na(match("errorlod",names(cross$geno[[i]])))) 
      stop("You first need to run calc.errorlod.")

    el <- cross$geno[[i]]$errorlod

    if(any(el > cutoff)) {
      o <- (el > cutoff)
      mar <- c(mar,colnames(el)[col(el)[o]])
      ind <- c(ind,row(el)[o])
      lod <- c(lod,el[o])
      chr <- c(chr,rep(names(cross$geno)[i],sum(o)))
      flag <- 1
    }
  }
  if(!flag) {
    if(msg) cat("\tNo errorlods above cutoff.\n")
    return(invisible(NULL))
  }
  o <- data.frame(chr=chr,ind=ind,marker=mar,errorlod=lod)[order(-lod,ind),]
  rownames(o) <- 1:nrow(o)
  o
}

# end of errorlod.R
