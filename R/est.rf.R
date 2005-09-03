######################################################################
#
# est.rf.R
#
# copyright (c) 2001-5, Karl W Broman, Johns Hopkins University
# last modified May, 2005
# first written Apr, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: est.rf, plot.rf, checkrf
#
######################################################################

######################################################################
#
# est.rf: Estimate sex-averaged recombination fractions between
#         all pairs of markers
#
######################################################################

est.rf <-
function(cross, maxit=4000, tol=1e-4) 
{
  n.chr <- nchr(cross)
  n.mar <- totmar(cross)
  n.ind <- nind(cross)
  mar.names <- unlist(lapply(cross$geno,function(a) colnames(a$data)))
  
  type <- class(cross)[1]
  chrtype <- sapply(cross$geno,class)

  xchrcol <- NULL
  fixX <- FALSE
  Geno <- NULL
  # create full genotype matrix
  for(i in 1:n.chr) {
    temp <- cross$geno[[i]]$data

    # treat X chromosome specially in an intercross
    if((type=="f2" || type=="f2ss") && chrtype[i]=="X") {
      fixX <- TRUE
      if(i != 1) xchrcol <- c(xchrcol,ncol(Geno)+(1:ncol(cross$geno[[i]]$data)))
      else xchrcol <- 1:ncol(cross$geno[[i]]$data)
      xchr <- temp
      xchr[is.na(xchr)] <- 0
      temp <- reviseXdata(type,"simple",getsex(cross),geno=temp)
    }
    Geno <- cbind(Geno,temp)
  }

  # which type of cross is this?
  if(type == "f2" || type=="f2ss") 
    cfunc <- "est_rf_f2"
  else if(type == "bc" || type=="risib" || type=="riself") 
    cfunc <- "est_rf_bc"
  else if(type == "4way") 
    cfunc <- "est_rf_4way"
  else {
    err <- paste("est.rf not available for cross type",
                 type, ".")
    stop(err)
  }

  Geno[is.na(Geno)] <- 0
  
  if(type=="bc" || type=="risib" || type=="riself")
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.mar),         # number of markers
            as.integer(Geno),
            rf = as.double(rep(0,n.mar*n.mar)),
            PACKAGE="qtl")
  else
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.mar),         # number of markers
            as.integer(Geno),
            rf = as.double(rep(0,n.mar*n.mar)),
            as.integer(maxit),
            as.double(tol),
            PACKAGE="qtl")

  cross$rf <- matrix(z$rf,ncol=n.mar)
  dimnames(cross$rf) <- list(mar.names,mar.names)

  if(fixX) {
    zz <- .C("est_rf_bc",
             as.integer(n.ind),
             as.integer(ncol(xchr)),
             as.integer(xchr),
             rf=as.double(rep(0,ncol(xchr)^2)),
             PACKAGE="qtl")
    zz <- matrix(zz$rf,ncol=ncol(xchr))
    cross$rf[xchrcol,xchrcol] <- zz
  }

  checkrf(cross, 3)
  cross
}

  

plot.rf <-
function(x, chr, which=c("both","lod","rf"), ...)
{
  which <- match.arg(which)
  
  if(!missing(chr)) x <- subset(x,chr=chr)
  
  if(is.na(match("rf",names(x)))) stop("You must run est.rf first.")
  g <- x$rf
  
  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=TRUE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))

  # if any of the rf's are NA (ie no data), put NAs in corresponding LODs
  if(any(is.na(g))) g[is.na(t(g))] <- NA

  # convert rf to -2*(log2(rf)+1); place 12's on the diagonal;
  #    anything above 12 replaced by 12;
  #    NA's replaced by -1
  g[row(g) > col(g) & g > 0.5] <- 0.5
  g[row(g) > col(g)] <- -4*(log2(g[row(g) > col(g)])+1)
  diag(g) <- 12
  g[!is.na(g) & g>12] <- 12
  
  g[is.na(g)] <- -1

  if(which=="lod") { # plot LOD scores 
    # copy upper triangle (LODs) to lower triangle (rec fracs)
    g[row(g) > col(g)] <- t(g)[row(g) > col(g)]
  }
  else if(which=="rf") { # plot recombination fractions
    # copy lower triangle (rec fracs) to upper triangle (LODs)
    g[row(g) < col(g)] <- t(g)[row(g) < col(g)]
  }
  br <- c(-1, seq(-1e-6, 12, length=65))


  image(1:ncol(g),1:nrow(g),t(g),ylab="Markers",xlab="Markers",breaks=br,
        col=c("lightgray",rev(rainbow(64,start=0,end=2/3))))
  
  # plot lines at the chromosome boundaries
  n.mar <- nmar(x)
  n.chr <- nchr(x)
  a <- c(0.5,cumsum(n.mar)+0.5)
  abline(v=a,xpd=FALSE)
  abline(h=a,xpd=FALSE)

  # this line adds a line above the image
  #     (the image function leaves it out)
  abline(h=0.5+nrow(g),xpd=FALSE)
  abline(v=0.5+nrow(g),xpd=FALSE)

  # add chromosome numbers
  a <- par("usr")
  wh <- cumsum(c(0.5,n.mar))
  for(i in 1:n.chr) 
    text(mean(wh[i+c(0,1)]),a[4]+(a[4]-a[3])*0.025,names(x$geno)[i])
  for(i in 1:n.chr) 
    text(a[2]+(a[2]-a[1])*0.025,mean(wh[i+c(0,1)]),names(x$geno)[i])

  if(which=="lod") title(main="Pairwise LOD scores")
  else if(which=="rf") title(main="Recombination fractions")
  else title("Pairwise recombination fractions and LOD scores")

  invisible()
}

######################################################################
# check for apparent errors in the recombination fractions
######################################################################
checkrf <-
function(cross, threshold=3)
{
  rf <- cross$rf
  n.mar <- nmar(cross)
  map <- pull.map(cross)
  n <- ncol(rf)
  mnam <- colnames(rf)
  whpos <- unlist(lapply(map,function(a) 1:length(a)))
  whchr <- rep(names(map),sapply(map,length))

  # first check whether a locus has "significant" pairwise recombination
  #     with rf > 0.5
  for(i in 1:n) {
    if(i == 1) {
      lod <- rf[1,-1]
      r <- rf[-1,1]
    }
    else if(i == n) {
      lod <- rf[-n,n]
      r <- rf[n,-n]
    }
    else {
      lod <- c(rf[1:(i-1),i],rf[i,(i+1):n])
      r <- c(rf[i,1:(i-1)],rf[(i+1):n,i])
    }

    # if rf > 1/2 and LOD > threshold for more than two other markers
    if(sum(!is.na(lod) & !is.na(r) & lod > threshold & r > 0.5) >= 2)
      warning("Genotypes potentially switched for marker ", mnam[i],
          paste(" (",whpos[i],")",sep=""), " on chr ", whchr[i], "\n")
    
  }

}

# end of est.rf.R
