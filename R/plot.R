######################################################################
#
# plot.R
#
# copyright (c) 2000-2001, Karl W Broman, Johns Hopkins University
# Oct, 2001; Sept, 2001; July, 2001; Apr, 2001; Feb, 2001; Mar, 2000
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: plot.missing, plot.map, plot.cross, plot.geno, plot.info
#
######################################################################

plot.missing <-
function(x,chr,reorder=FALSE,main="Missing genotypes",...) 
{
  cross <- x
  if(!missing(chr)) cross <- pull.chr(cross,chr)
  
  # get full genotype data into one matrix
  Geno <- cross$geno[[1]]$data
  if(length(cross$geno) > 1) 
    for(i in 2:length(cross$geno))
      Geno <- cbind(Geno,cross$geno[[i]]$data)

  # reorder the individuals according to their phenotype
  o <- 1:nrow(Geno)
  if(reorder) {
    # if reorder is a number, use the corresponding phenotype
    if(is.numeric(reorder)) o <- order(cross$pheno[,reorder])

    # otherwise, order according to the sum of the phenotypes
    else o <- order(apply(cross$pheno,1,sum))
  }

  # make matrix with  0 where genotype data is missing
  #                   1 where data is not missing
  #                 0.5 where data is partially missing
  type <- class(cross)[1]
  g <- t(Geno[o,])
  g[is.na(g)] <- 0
  if(type == "bc") {
    g[g > 0] <- 1
  }
  else if(type=="f2") {
    g[g > 0 & g < 4] <- 1
    g[g > 3] <- 0.5
  }
  else if(type=="4way") {
    g[g > 0 & g < 5] <- 1
    g[g > 4] <- 0.5
  }
  else {
    g[g > 0] <- 1
  }

  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=TRUE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))

  colors <- c("#000000", "gray80", "#FFFFFF")

  # plot grid with black pixels where there is missing data
  image(1:nrow(g),1:ncol(g),g,ylab="Individuals",xlab="Markers",col=colors,zlim=c(0,1))

  # plot lines at the chromosome boundaries
  n.mar <- nmar(cross)
  a <- c(0.5,cumsum(n.mar)+0.5)
#  abline(v=a)
  # the following makes the lines go slightly above the plotting region
  b <- par("usr")
  segments(a,b[3],a,b[4]+diff(b[3:4])*0.02)

  # this line adds a line above the image
  #     (the image function seems to leave it out)
  abline(h=0.5+c(0,ncol(g)),xpd=FALSE)

  # add chromosome numbers
  a <- par("usr")
  wh <- cumsum(c(0.5,n.mar))
  for(i in 1:length(n.mar)) 
    text(mean(wh[i+c(0,1)]),a[4]+(a[4]-a[3])*0.025,names(cross$geno)[i])

  title(main=main)
}

plot.map <-
function(x,map2,horizontal=FALSE,...) 
{
  map <- x
  # figure out if the input is a cross (containing a map)
  #    or is the map itself
  if(!is.na(match("geno",names(map)))) 
    map <- pull.map(map)

  sex.sp <- FALSE

  if(is.matrix(map[[1]])) { # sex-specific maps
    one.map <- FALSE
    sex.sp <- TRUE
    if(!missing(map2)) {
      if(is.logical(map2)) {
        horizontal <- map2
        map2 <- lapply(map,function(a) a[2,])
        map <- lapply(map,function(a) a[1,])
      }
      else {
        if(!is.na(match("geno",names(map2))))
          map2 <- pull.map(map2)
        Map1 <- lapply(map,function(a) a[1,])
        Map2 <- lapply(map,function(a) a[2,])
        Map3 <- lapply(map2,function(a) a[1,])
        Map4 <- lapply(map2,function(a) a[2,])
        old.mfrow <- par("mfrow")
        on.exit(par(mfrow=old.mfrow))
        par(mfrow=c(2,1))
        plot.map(Map1,Map3,horizontal)
        plot.map(Map2,Map4,horizontal)
        return(invisible())
      }
    }
    else {
      map2 <- lapply(map,function(a) a[2,])
      map <- lapply(map,function(a) a[1,])
    }
  }
  else {
    # determine whether a second map was given
    if(!missing(map2)) {
      if(is.logical(map2)) { # assume "map2" should be "horizontal"
        horizontal <- map2
        map2 <- NULL
        one.map <- TRUE
      }
      else { # determine if it is a cross object
        if(!is.na(match("geno",names(map2))))
          map2 <- pull.map(map2)
        one.map <- FALSE
      }
    }
    else one.map <- TRUE
  }
       
  if(one.map) {
    n.chr <- length(map)
    map <- lapply(map, function(a) a-min(a))
    maxlen <- max(unlist(lapply(map,max)))

    if(horizontal) {
      old.xpd <- par("xpd")
      old.yaxt <- par("yaxt")
      par(xpd=TRUE,yaxt="n")
      on.exit(par(xpd=old.xpd,yaxt=old.yaxt))
      
      plot(0,0,type="n",xlim=c(0,maxlen),ylim=c(0.5,n.chr+0.5),
	   xlab="Location (cM)", ylab="Chromosome")
      a <- par("usr")
      
      for(i in 1:n.chr) {
	lines(c(min(map[[i]]),max(map[[i]])),n.chr+1-c(i,i))
	nmar <- length(map[[i]])
	for(j in 1:nmar)
	  lines(rep(map[[i]][j],2),n.chr+1-i+c(-1/4,1/4))

	# add chromosome label
	text(a[1]-(a[2]-a[1])*0.02,n.chr+1-i,names(map)[i],adj=1)
	lines(c(a[1],a[1]-(a[2]-a[1])*0.01),rep(n.chr+1-i,2))
      }
    }
    else {
      old.xpd <- par("xpd")
      old.xaxt <- par("xaxt")
      old.las <- par("las")
      par(xpd=TRUE,xaxt="n",las=1)
      on.exit(par(xpd=old.xpd,xaxt=old.xaxt,las=old.las))
      
      plot(0,0,type="n",ylim=c(0,maxlen),xlim=c(0.5,n.chr+0.5),
	   ylab="Location (cM)", xlab="Chromosome")
      
      a <- par("usr")
      
      for(i in 1:n.chr) {
	lines(c(i,i), c(min(map[[i]]),max(map[[i]])))
	nmar <- length(map[[i]])
	for(j in 1:nmar)
	  lines(i+c(-1/4,1/4),rep(map[[i]][j],2))

        # add chromosome label
	text(i,a[3]-(a[4]-a[3])*0.04,names(map)[i])
	lines(rep(i,2),c(a[3],a[3]-(a[4]-a[3])*0.02))
      }
    }
    title(main="Genetic map")
  }
  else {
    # check that maps conform
    if(is.matrix(map2[[1]]))
      stop("Second map appears to be a sex-specific map.")
    if(length(map) != length(map2))
      stop("Maps have different numbers of chromosomes.")
    if(any(sapply(map,length) != sapply(map2,length)))
      stop("Maps have different numbers of markers.")

    map1 <- lapply(map,function(a) a-a[1])
    map2 <- lapply(map2,function(a) a-a[1])

    n.chr <- length(map1)
    maxloc <- max(c(unlist(lapply(map1,max)),unlist(lapply(map2,max))))

    if(!horizontal) {
      old.xpd <- par("xpd")
      old.xaxt <- par("xaxt")
      old.las <- par("las")
      par(xpd=TRUE,xaxt="n",las=1)
      on.exit(par(xpd=old.xpd,xaxt=old.xaxt,las=old.las))

      plot(0,0,type="n",ylim=c(0,maxloc),xlim=c(0.5,n.chr+0.5),
           ylab="Location (cM)", xlab="Chromosome")

      a <- par("usr")
    
      for(i in 1:n.chr) {
      
        if(max(map2[[i]]) < max(map1[[i]])) 
          map2[[i]] <- map2[[i]] + (max(map1[[i]])-max(map2[[i]]))/2
        else 
          map1[[i]] <- map1[[i]] + (max(map2[[i]])-max(map1[[i]]))/2
        
        lines(c(i-0.3,i-0.3), c(min(map1[[i]]),max(map1[[i]])))
        lines(c(i+0.3,i+0.3), c(min(map2[[i]]),max(map2[[i]])))
        
        nmar <- length(map1[[i]])
        for(j in 1:nmar)
          lines(c(i-0.3,i+0.3),c(map1[[i]][j],map2[[i]][j]))

        # add chromosome label
        text(i,a[3]-(a[4]-a[3])*0.04,names(map1)[i])
        lines(rep(i,2),c(a[3],a[3]-(a[4]-a[3])*0.02))
      }
    }
    else {
      old.xpd <- par("xpd")
      old.yaxt <- par("yaxt")
      old.las <- par("las")
      par(xpd=TRUE,yaxt="n",las=1)
      on.exit(par(xpd=old.xpd,yaxt=old.yaxt,las=old.las))

      plot(0,0,type="n",xlim=c(0,maxloc),ylim=c(0.5,n.chr+0.5),
           xlab="Location (cM)", ylab="Chromosome")

      a <- par("usr")
    
      for(i in 1:n.chr) {
      
        if(max(map2[[i]]) < max(map1[[i]])) 
          map2[[i]] <- map2[[i]] + (max(map1[[i]])-max(map2[[i]]))/2
        else 
          map1[[i]] <- map1[[i]] + (max(map2[[i]])-max(map1[[i]]))/2
        
        lines(c(min(map2[[i]]),max(map2[[i]])), c(n.chr-i-0.3+1,n.chr-i+1-0.3))
        lines(c(min(map1[[i]]),max(map1[[i]])), c(n.chr-i+1+0.3,n.chr+1-i+0.3))
        
        nmar <- length(map1[[i]])
        for(j in 1:nmar)
          lines(c(map2[[i]][j],map1[[i]][j]), c(n.chr+1-i-0.3,n.chr+1-i+0.3))

        # add chromosome label
        text(a[1]-diff(a[1:2])*0.04,n.chr+1-i, names(map1)[i])
        lines(c(a[1],a[1]-diff(a[1:2])*0.02), rep(n.chr+1-i,2))
      }

    }
    if(!sex.sp) title(main="Comparison of genetic maps")
    else title(main="Genetic map")
  }    

}


plot.cross <-
function(x,...)
{
  cross <- x
  
  old.yaxt <- par("yaxt")
  if(ncol(cross$pheno) > 2) {
    old.ask <- par("ask")
    par(ask=TRUE)
    on.exit(par(ask=old.ask,yaxt=old.yaxt))
  }
  else {
    old.mfrow <- par("mfrow")
    par(mfrow=c(2,2))
    on.exit(par(mfrow=old.mfrow,yaxt=old.yaxt))
  }

  plot.missing(cross)
  plot.map(cross)
  par(yaxt = "n")
  for(i in 1:ncol(cross$pheno)) 
    hist(cross$pheno[,i],nclass=round(sqrt(nrow(cross$pheno))+5),
	 xlab=colnames(cross$pheno)[i],prob=TRUE, ylab="",
	 main=paste("Histogram of", colnames(cross$pheno)[i]))
}


##################################################r####################
#
# plot.geno: Plot genotypes for a specified chromosome, with likely
#           genotyping errors indicated. 
#
######################################################################

plot.geno <-
function(x, chr, ind, horizontal=FALSE, cutoff=2,
         method=c("lod","argmax"),min.sep=1,...)
{
  cross <- x  
  method <- match.arg(method)
  cross <- pull.chr(cross,chr)
  type <- class(cross)[1]
  
  if(type != "bc" && type != "f2")
    stop("This function has only been coded for bc and f2 crosses.")

  if(method=="lod") {
    if(is.na(match("errorlod",names(cross$geno[[1]])))) {
      warning("First running calc.errorlod.")
      cross <- calc.errorlod(cross,error.prob=0.01)
    }
  }
  else {
    if(is.na(match("argmax",names(cross$geno[[1]])))) {
      warning("First running argmax.geno.")
      cross <- argmax.geno(cross,error.prob=0.01)
    }
  }
  
  # if necessary, discard parts of argmax that are not at markers
  if(method=="argmax") {
    wh <- grep("^loc\-*[0-9]+",colnames(cross$geno[[1]]$argmax))
    if(length(wh) > 0) 
      cross$geno[[1]]$argmax <- cross$geno[[1]]$argmax[,-wh]
  }

  # indicators for apparent errors
  errors <- matrix(0,ncol=ncol(cross$geno[[1]]$data),
                   nrow=nrow(cross$geno[[1]]$data))
  dimnames(errors) <- dimnames(cross$geno[[1]]$data)

  if(method=="lod") top <- top.errorlod(cross,1,cutoff,FALSE)
  else top <- find.errors(cross,1,msg=FALSE)
  if(length(top) > 0)
    for(i in 1:nrow(top))
      errors[top[i,2],as.character(top[i,3])] <- 1

  # map, data, errors
  map <- cross$geno[[1]]$map
  if(is.matrix(map)) map <- map[1,] # if sex-specific map
  L <- diff(range(map))
  min.d <- L*min.sep/100
  d <- diff(map)
  d[d < min.d] <- min.d
  map <- cumsum(c(0,d))

  data <- cross$geno[[1]]$data
  if(!missing(ind)) {
    data <- data[ind,]
    errors <- errors[ind,]
  }
  n.ind <- nrow(errors)

  color <- c("white","gray60","black","green","orange","red")

  if(horizontal==TRUE) {
    plot(0,0,type="n",xlab="Position (cM)",ylab="Individual",
         main=paste("Chromosome",names(cross$geno)[1]),
         ylim=c(0.5,n.ind+0.5),xlim=c(0,max(map)))
    segments(0,1:n.ind,max(map),1:n.ind)

    # AA genotypes
    tind <- rep(1:n.ind,length(map));tind[is.na(data)] <- NA
    ind <- tind; ind[!is.na(data) & data!=1] <- NA
    x <- rep(map,rep(n.ind,length(map)))
    points(x,ind,pch=16,col=color[1])
    points(x,ind,pch=1)

    # AB genotypes
    ind <- tind; ind[!is.na(data) & data!=2] <- NA
    if(type=="bc")
      points(x,ind,pch=16,col=color[3]) 
    else {
      points(x,ind,pch=16,col=color[2])
      points(x,ind,pch=1)
    }

    if(type=="f2") {
      # BB genotypes
      ind <- tind; ind[!is.na(data) & data!=3] <- NA
      points(x,ind,pch=16,col=color[3])

      # not BB (D in mapmaker/qtl) genotypes
      ind <- tind; ind[!is.na(data) & data!=4] <- NA
      points(x,ind,pch=16,col=color[4])
      points(x,ind,pch=1)

      # not AA (C in mapmaker/qtl) genotypes
      ind <- tind; ind[!is.na(data) & data!=5] <- NA
      points(x,ind,pch=16,col=color[5])
      points(x,ind,pch=1)
    }

    # plot map
    u <- par("usr")
    segments(map,u[3],map,(u[3]+1)/2)
    segments(map,u[4],map,(n.ind+u[4])/2)

    if(any(errors)) {
      ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
      points(x,ind,pch=0,col=color[6],cex=1.5)
    }

  }
  else {
    plot(0,0,type="n",ylab="Position (cM)",xlab="Individual",
         main=paste("Chromosome",names(cross$geno)[1]),
         xlim=c(0.5,n.ind+0.5),ylim=c(0,max(map)))
    segments(1:n.ind,0,1:n.ind,max(map))
    
    # AA genotypes
    tind <- rep(1:n.ind,length(map));tind[is.na(data)] <- NA
    ind <- tind; ind[!is.na(data) & data!=1] <- NA
    y <- rep(map,rep(n.ind,length(map)))
    points(ind,y,pch=16,col="white")
    points(ind,y,pch=1)

    # AB genotypes
    ind <- tind; ind[!is.na(data) & data!=2] <- NA
    if(type=="bc")
      points(ind,y,pch=16,col=color[3])
    else {
      points(ind,y,pch=16,col=color[2])
      points(ind,y,pch=1)
    }

    if(type=="f2") {
      # BB genotypes
      ind <- tind; ind[!is.na(data) & data!=3] <- NA
      points(ind,y,pch=16,col=color[3])

      # not BB genotypes
      ind <- tind; ind[!is.na(data) & data!=4] <- NA
      points(ind,y,pch=16,col=color[4])
      points(ind,y,pch=1)

      # not AA genotypes
      ind <- tind; ind[!is.na(data) & data!=5] <- NA
      points(ind,y,pch=16,col=color[5])
      points(ind,y,pch=1)
    }

    # plot map
    u <- par("usr")
    segments(u[1],map,(u[1]+1)/2,map)
    segments(u[2],map,(n.ind+u[2])/2,map)

    if(any(errors)) {
      ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
      points(ind,y,pch=0,col=color[6],cex=1.5)
    }
  }
}
    
######################################################################
#
# plot.info: Plot the proportion of missing information in the
#            genotype data.
#
######################################################################

plot.info <-
function(x,chr,which=c("both","entropy","variance"),return.result=FALSE,...)
{
  cross <- x
  which <- match(match.arg(which),c("entropy","variance","both"))-1

  if(!missing(chr)) cross <- pull.chr(cross,chr)
  results <- NULL

  n.chr <- nchr(cross)
  if(is.na(match("prob",names(cross$geno[[1]])))) { # need to run calc.genoprob
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }
  gap <- attr(cross$geno[[1]]$prob,"off.end")*2+10 # gap between chr in plot
  n.ind <- nind(cross)
  for(i in 1:n.chr) {
    n.gen <- dim(cross$geno[[i]]$prob)[3]
    n.pos <- ncol(cross$geno[[i]]$prob)

    # calculate information (between 0 and 1)
    info <- .C("R_info",
               as.integer(n.ind),
               as.integer(n.pos),
               as.integer(n.gen),
               as.double(cross$geno[[i]]$prob),
               info1=as.double(rep(0,n.pos)),
               info2=as.double(rep(0,n.pos)),
               as.integer(which))

    if(which != 1) { # rescale entropy version
      if(n.gen==3) maxent <- 1.5*log(2)
      else maxent <- log(n.gen)
      info$info1 <- -info$info1/maxent
    }
    if(which != 0) { # rescale variance version
      maxvar <- c(0.25,0.5,1.25)[n.gen-1]
      info$info2 <- info$info2/maxvar
    }

    # reconstruct map
    map <- create.map(cross$geno[[i]]$map,
                      attr(cross$geno[[i]]$prob,"step"),
                      attr(cross$geno[[i]]$prob,"off.end"))
    if(is.matrix(map)) map <- map[1,]

    z <- data.frame(chr=rep(names(cross$geno)[i],length(map)),pos=map,
                    "Missing information"=info$info1,
                    "Missing information"=info$info2)
    w <- names(map)
    o <- grep("^loc\-*[0-9]+",w)
    if(length(o) > 0) 
      w[o] <- paste(w[o],names(cross$geno)[i],sep=".c")
    rownames(z) <- w
    results <- rbind(results,z)
  }

  if(which==0) plot.scanone(results,ylim=c(0,1),gap=gap)
  else if(which==1) plot.scanone(results[,-3],ylim=c(0,1),gap=gap)
  else if(which==2) plot.scanone(results,results[,-3],ylim=c(0,1),gap=gap)
  
  if(return.result) {
    colnames(results)[3:4] <- c("misinfo.entropy","misinfo.variance")
    class(results) <- c("scanone","data.frame")
    return(results)
  }
  else invisible()
}

# end of plot.R
