######################################################################
#
# est.rf.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# last modified Nov, 2001
# first written Apr, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: est.rf, plot.rf
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
  
  Geno <- NULL
  # create full genotype matrix
  for(i in 1:n.chr) 
    Geno <- cbind(Geno,cross$geno[[i]]$data)

  # which type of cross is this?
  type <- class(cross)[1]
  if(type == "f2") 
    cfunc <- "est_rf_f2"
  else if(type == "bc") 
    cfunc <- "est_rf_bc"
  else if(type == "4way") 
    cfunc <- "est_rf_4way"
  else stop(paste("est.rf not available for cross type",
                  type, "."))

  Geno[is.na(Geno)] <- 0
  
  if(type=="bc")
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
  
}

# end of est.rf.R
