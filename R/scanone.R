######################################################################
#
# scanone.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# Oct, 2001; Sept, 2001; May, 2001, Apr, 2001; Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: scanone, plot.scanone, scanone.perm
#           summary.scanone, print.summary.scanone
#
######################################################################

######################################################################
#
# scanone: scan genome, calculating LOD scores with single QTL model
#          (currently covariates not allowed; will be added later)
#
######################################################################

scanone <-
function(cross, chr, pheno.col=1, method=c("im","anova","hk","imp"),
         start=NULL, maxit=1000, tol=1e-8)
{
  method <- match.arg(method)
  if(method=="imp")
    warning("We don't have the imputation method working yet.")

  if(!missing(chr)) cross <- pull.chr(cross,chr)

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

  if(method == "im") {
    m <- mean(pheno)
    s <- sd(pheno)*sqrt((n.ind-1)/n.ind)
    llik0 <- sum(log10(dnorm(pheno,m,s)))
  }
  else if(method=="hk")
    lrss0 <- log10(sum((pheno-mean(pheno))^2))

  results <- NULL

  if(method=="im") {
    if(is.null(start)) std.start <- 1
    else if(length(start)==1) std.start <- -1
    else std.start <- 0
  }

  # calculate genotype probabilities one chromosome at a time
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
    else stop(paste("scanone not available for cross type",
                    type, "."))

    # starting values for interval mapping
    if(method=="im") {
      this.start <- rep(0,n.gen+1)
      if(std.start == 0) {
        if(length(start) < n.gen+1)
          stop(paste("Length of start argument should be 0, 1 or", n.gen+1))
        this.start <- c(start[1:n.gen],start[length(start)])
      }
    }

    # pull out reconstructed genotypes (anova)
    # or genotype probabilities (im or hk)

    if(method == "anova") {
      cfunc <- "R_scanone_anova"
      newgeno <- cross$geno[[i]]$data
      newgeno <- newgeno[keep.ind,]
      newgeno[is.na(newgeno)] <- 0 

      # discard partially informative genotypes
      if(type=="f2" || type=="f2ss") newgeno[newgeno>3] <- 0
      if(type=="4way") newgeno[newgeno>4] <- 0

      n.pos <- ncol(newgeno)
      map <- cross$geno[[i]]$map
      if(is.matrix(map)) map <- map[1,]
    }
    else if(method == "imp") {
      if(is.na(match("draws",names(cross$geno[[i]])))) { # need to run sim.geno
        warning("First running sim.geno.")
        cross <- sim.geno(cross)
      }

      draws <- cross$geno[[i]]$draws
      n.pos <- ncol(draws)
      n.draws <- dim(draws)[3]
      draws <- draws[keep.ind,,]

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$draws,"step"),
                        attr(cross$geno[[i]]$draws,"off.end"))
      if(is.matrix(map)) map <- map[1,]

      cfunc <- "R_scanone_imp"
    }
    else {
      if(is.na(match("prob",names(cross$geno[[i]])))) { # need to run calc.genoprob
        warning("First running calc.genoprob.")
        cross <- calc.genoprob(cross)
      }
      genoprob <- cross$geno[[i]]$prob
      n.pos <- ncol(genoprob)
      genoprob <- genoprob[keep.ind,,]

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$prob,"step"),
                        attr(cross$geno[[i]]$prob,"off.end"))
      if(is.matrix(map)) map <- map[1,]

      if(method == "im") cfunc <- "R_scanone_im"
      else cfunc <- "R_scanone_hk"
    }

    # call the C function
    if(method == "anova") {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.integer(newgeno),       # genotype data
              as.double(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+2))),
              PACKAGE="qtl")
    }
    else if(method=="imp") {
      z <- .C(cfunc,
              as.integer(n.ind),
              as.integer(n.pos),
              as.integer(n.gen),
              as.integer(n.draws),
              as.integer(draws),
              as.double(pheno),
              result=as.double(rep(0,n.pos)),
              PACKAGE="qtl")
    }
    else if(method=="hk") {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+2))),
              PACKAGE="qtl")
    }
    else { # interval mapping
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(pheno),          # phenotype data
              result=as.double(rep(0,n.pos*(n.gen+2))),
              as.integer(std.start),
              as.double(this.start),
              as.integer(maxit),
              as.double(tol),
              PACKAGE="qtl")
    }

    z <- matrix(z$result,nrow=n.pos)

    if(method == "im")
      z[,1] <- z[,1] - llik0
    else if(method == "hk")
      z[,1] <- (n.ind/2)*(lrss0-z[,1])
#    else if(method == "imp")
#      z[,1] <- z[,1] - (n.gen-1)/2*log10(n.ind)

    if(method != "imp") {
      if(type=="f2" && class(cross$geno[[i]])=="X") # add BB column
        z <- cbind(z[,1:3],rep(NA,n.pos),z[,4])

      colnames(z) <- c("lod",gen.names,"sigma")
    }
    else colnames(z) <- c("lod")
      
    w <- names(map)
    o <- grep("^loc\-*[0-9]+",w)
    if(length(o) > 0) 
      w[o] <- paste(w[o],names(cross$geno)[i],sep=".c")
    rownames(z) <- w
    
    z <- as.data.frame(z)
    z <- cbind(chr=rep(names(cross$geno)[i],length(map)), pos=map, z)
    rownames(z) <- w
    results <- rbind(results,z)
  }

  # replace any lod = NaN with 0
  results[is.na(results[,3]),3] <- 0

  class(results) <- c("scanone",class(results))
  results
}

  

######################################################################
#
# plot.scanone: plot output from scanone
#
######################################################################

plot.scanone <- 
function(x,output2,output3,chr,incl.markers=TRUE,ylim,
         lty=c(1,2,3),col="black",lwd=2,add=FALSE,gap=25,...)
{
  output <- x
  second <- third <- TRUE
  if(missing(output2) && missing(output3)) 
     second <- third <- FALSE
  if(missing(output3))
    third <- FALSE
  if(missing(output2))
    second <- FALSE

  if(length(lty)==1) lty <- rep(lty,3)
  if(length(lwd)==1) lwd <- rep(lwd,3)
  if(length(col)==1) col <- rep(col,3)

  # pull out desired chromosomes
  if(missing(chr))
    chr <- unique(as.character(output[,1]))

  if(length(chr) == 0) chr <- sort(unique(output[,1]))
  else if(all(chr < 0)) { 
    a <- sort(unique(output[,1]))
    chr <- a[-match(-chr,a)]
  }
  output <- output[!is.na(match(output[,1],chr)),]
  if(second) output2 <- output2[!is.na(match(output2[,1],chr)),]
  if(third) output3 <- output3[!is.na(match(output3[,1],chr)),]
  
  # beginning and end of chromosomes
  temp <- grep("^loc\-*[0-9]+",rownames(output))
  if(length(temp)==0) temp <- output
  else temp <- output[-temp,]
  begend <- matrix(unlist(tapply(temp[,2],temp[,1],range)),ncol=2,byrow=TRUE)
  len <- begend[,2]-begend[,1]

  # locations to plot start of each chromosome
  start <- gap/2+c(0,cumsum(len+gap))-c(begend[,1],0)

  maxx <- sum(len+gap)
  maxy <- max(output[,3])
  if(second) maxy <- max(c(maxy,output2[,3]))
  if(third) maxy <- max(c(maxy,output3[,3]))

  # graphics parameters
  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=TRUE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))

  # make frame of plot
  if(missing(ylim)) ylim <- c(0,maxy)

  if(!add)
    plot(0,0,ylim=ylim,xlim=c(0,maxx),type="n",
         xlab="Map position (cM)",ylab=dimnames(output)[[2]][3])

  for(i in 1:length(chr)) {
    # plot first output
    x <- output[output[,1]==chr[i],2]+start[i]
    y <- output[output[,1]==chr[i],3]
    lines(x,y,lwd=lwd[1],lty=lty[1],col=col[1])

    # plot chromosome number
    a <- par("usr")
    if(!add) {
      tloc <- mean(c(start[i],start[i+1]-gap))
      text(tloc,a[4]+(a[4]-a[3])*0.03,as.character(chr[i]))
      lines(rep(tloc,2),c(a[4],a[4]+(a[4]-a[3])*0.015))
    }

    # plot second output
    if(second) {
      x <- output2[output2[,1]==chr[i],2]+start[i]
      y <- output2[output2[,1]==chr[i],3]
      lines(x,y,lty=lty[2],col=col[2],lwd=lwd[2])
    }

    if(third) {
      x <- output3[output3[,1]==chr[i],2]+start[i]
      y <- output3[output3[,1]==chr[i],3]
      lines(x,y,lty=lty[3],col=col[3],lwd=lwd[3])
    }

    # plot lines at marker positions
    if(incl.markers && !add) {
      nam <- dimnames(output)[[1]][output[,1]==chr[i]]
      wh.genoprob <- (1:length(nam))[grep("^loc\-*[0-9]+",nam)]
      if(length(wh.genoprob)==0) wh.genoprob <- 1:length(nam)
      else wh.genoprob <- (1:length(nam))[-wh.genoprob]
      pos <- output[output[,1]==chr[i],2][wh.genoprob]+start[i]
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
function(cross, chr, pheno.col=1, method=c("im","anova","hk"),
         start=NULL,n.perm=1000, maxit=1000, tol=1e-8)
{
  method <- match.arg(method)

  if(method=="hk") {
    warning("We don't have Haley-Knott working yet; running IM instead.")
    method <- "im"
  }

  if(!missing(chr))
    cross <- pull.chr(cross,chr)

  # check phenotypes
  if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
  if(pheno.col < 1 || pheno.col > nphe(cross))
    stop("Specified phenotype column is invalid.")

  if(method=="im") {
    if(is.null(start)) std.start <- 1
    else if(length(start)==1) std.start <- -1
    else std.start <- 0
  }

  pheno <- cross$pheno[,pheno.col]
  keep.ind <- (1:length(pheno))[!is.na(pheno)]
  pheno <- pheno[keep.ind]
  n.ind <- length(keep.ind)
  n.chr <- nchr(cross)
  type <- class(cross)[1]

  if(method == "im") {
    m <- mean(pheno)
    s <- sd(pheno)*sqrt((n.ind-1)/n.ind)
    llik0 <- sum(log10(dnorm(pheno,m,s)))
  }
  else lrss0 <- log10(sum((pheno-mean(pheno))^2))

  results <- NULL

  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {

    # which type of cross is this?
    if(type == "f2") {
      if(class(cross$geno[[i]]) == "A")  # autosomal
        n.gen <- 3
      else                              # X chromsome 
        n.gen <- 2
    }
    else if(type == "bc") 
      n.gen <- 2
    else if(type == "4way") 
      n.gen <- 4
    else stop(paste("scanone not available for cross type",
                    type, "."))

    if(method=="im") {
      # starting values
      this.start <- rep(0,n.gen+1)
      if(std.start == 0) {
        if(length(start) < n.gen+1)
          stop(paste("Length of start argument should be 0, 1 or", n.gen+1))
        this.start <- c(start[1:n.gen],start[length(start)])
      }
    }

    # pull out genotypes (anova) or genotype probabilities (im or hk)

    if(method == "anova") {
      cfunc <- "scanone_anova_perm"
      newgeno <- cross$geno[[i]]$data
      newgeno <- newgeno[keep.ind,]
      newgeno[is.na(newgeno)] <- 0 

      # discard partially informative genotypes
      if(type=="f2" || type=="f2ss") newgeno[newgeno>3] <- 0
      if(type=="4way") newgeno[newgeno>4] <- 0

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

      map <- create.map(cross$geno[[i]]$map,
                        attr(cross$geno[[i]]$prob,"step"),
                        attr(cross$geno[[i]]$prob,"off.end"))
      if(is.matrix(map)) map <- map[1,]

      if(method == "im") cfunc <- "scanone_im_perm"
      else cfunc <- "scanone_hk_perm"
    }

    # call the C function
    if(method=="anova") 
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.integer(newgeno),       # genotype data
              as.double(pheno),          # phenotype data
              as.integer(n.perm),
              result=as.double(rep(0,n.perm)),
              PACKAGE="qtl")
    else if(method=="hk") 
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(pheno),          # phenotype data
              as.integer(n.perm),
              result=as.double(rep(0,n.perm)),
              PACKAGE="qtl")
    else # interval mapping
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(n.gen),         # number of possible genotypes
              as.double(genoprob),       # genotype probabilities
              as.double(pheno),          # phenotype data
              as.integer(n.perm),
              result=as.double(rep(0,n.perm)),
              as.integer(std.start),
              as.double(this.start),
              as.integer(maxit),
              as.double(tol),
              PACKAGE="qtl")

    if(method == "im")
      z <- z$result - llik0
    else 
      z <- (n.ind/2)*(lrss0-z$result)


    results <- cbind(results,z)
  }

  results <- cbind(apply(results,1,max),results)
  colnames(results) <- c("max",names(cross$geno))

  results
}


# give, for each chromosome, the output line at the maximum LOD
summary.scanone <-
function(object,threshold=0,...)
{
  output <- lapply(split(object,object[,1]),
                   function(b) b[b[,3]==max(b[,3]),])
  results <- output[[1]]
  if(length(output) > 1)
    for(i in 2:length(output))
      results <- rbind(results,output[[i]])
  class(results) <- c("summary.scanone","data.frame")
  if(!any(results[,3] >= threshold)) {
    cat("    There were no LOD peaks above the threshold.\n")
    invisible()
  }
  else {
    return(results[results[,3] >= threshold,])
  }
}

# print output of summary.scanone
print.summary.scanone <-
function(x,...)
{
  x[,-(1:2)] <- round(data.frame(x[,-(1:2)]),6)
  print.data.frame(x,digits=2)
}


# end of scanone.R
