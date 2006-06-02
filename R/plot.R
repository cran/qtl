######################################################################
#
# plot.R
#
# copyright (c) 2000-6, Karl W Broman, Johns Hopkins University
#       [modifications of plot.cross from Brian Yandell]
# last modified Jun, 2006
# first written Mar, 2000
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: plot.missing, plot.map, plot.cross, plot.geno, plot.info,
#           plot.pxg
#
######################################################################

plot.missing <-
function(x, chr, reorder=FALSE, main="Missing genotypes", ...) 
{
  cross <- x
  if(!missing(chr)) cross <- subset(cross,chr=chr)
  
  # get full genotype data into one matrix
  Geno <- cross$geno[[1]]$data
  if(length(cross$geno) > 1) 
    for(i in 2:length(cross$geno))
      Geno <- cbind(Geno,cross$geno[[i]]$data)

  # reorder the individuals according to their phenotype
  o <- 1:nrow(Geno)
  if(reorder) {
    # if reorder is a number, use the corresponding phenotype
    if(is.numeric(reorder)) {
      if(reorder < 1 || reorder > nphe(cross)) 
        stop("reorder should be TRUE, FALSE, or an integer between 1 and", nphe(cross))

      o <- order(cross$pheno[,reorder])
    }

    # otherwise, order according to the sum of the phenotypes
    else o <- order(apply(cross$pheno,1,sum))
  }

  # make matrix with  0 where genotype data is missing
  #                   1 where data is not missing
  #                 0.5 where data is partially missing
  type <- class(cross)[1]
  g <- t(Geno[o,])
  g[is.na(g)] <- 0
  if(type == "bc" || type=="risib" || type=="riself") 
    g[g > 0] <- 1
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
  n.chr <- nchr(cross)
  a <- c(0.5,cumsum(n.mar)+0.5)

  # the following makes the lines go slightly above the plotting region
  b <- par("usr")
  segments(a,b[3],a,b[4]+diff(b[3:4])*0.02)

  # this line adds a line above the image
  #     (the image function seems to leave it out)
  abline(h=0.5+c(0,ncol(g)),xpd=FALSE)

  # add chromosome numbers
  a <- par("usr")
  wh <- cumsum(c(0.5,n.mar))
  for(i in 1:n.chr)
    text(mean(wh[i+c(0,1)]),a[4]+(a[4]-a[3])*0.025,names(cross$geno)[i])

  title(main=main)
  invisible()
}

plot.map <-
function(x, map2, chr, horizontal=FALSE, shift=TRUE, ...) 
{
  map <- x
  # figure out if the input is a cross (containing a map)
  #    or is the map itself
  if(length(class(map))>1 && class(map)[2]=="cross")
    map <- pull.map(map)
  if(!missing(map2) && length(class(map2))>1 && class(map2)[2]=="cross")
    map2 <- pull.map(map2)
  
  if(!missing(chr)) {
    map <- map[chr]
    if(!missing(map2)) map2 <- map2[chr]
  }

  sex.sp <- FALSE

  if(is.matrix(map[[1]])) { # sex-specific map
    one.map <- FALSE
    sex.sp <- TRUE
    if(!missing(map2)) {
      if(is.logical(map2)) {
        horizontal <- map2
        map2 <- lapply(map,function(a) a[2,])
        map <- lapply(map,function(a) a[1,])
      }
      else {
        Map1 <- lapply(map,function(a) a[1,,drop=TRUE])
        Map2 <- lapply(map,function(a) a[2,,drop=TRUE])
        Map3 <- lapply(map2,function(a) a[1,,drop=TRUE])
        Map4 <- lapply(map2,function(a) a[2,,drop=TRUE])
        old.mfrow <- par("mfrow")
        on.exit(par(mfrow=old.mfrow))

        par(mfrow=c(2,1))
        plot.map(Map1,Map3,horizontal=horizontal,shift=shift)
        plot.map(Map2,Map4,horizontal=horizontal,shift=shift)
        return(invisible())
      }
    }
    else {
      map2 <- lapply(map,function(a) a[2,])
      map <- lapply(map,function(a) a[1,])
    }
  }
  else { # single map
    # determine whether a second map was given
    if(!missing(map2)) 
      one.map <- FALSE
    else one.map <- TRUE
  }
       
  if(one.map) {
    n.chr <- length(map)
    if(shift) map <- lapply(map, function(a) a-a[1])
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
      
      plot(0,0,type="n",ylim=c(maxlen,0),xlim=c(0.5,n.chr+0.5),
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

    map1 <- map
    if(shift) {
      map1 <- lapply(map1,function(a) a-a[1])
      map2 <- lapply(map2,function(a) a-a[1])
    }

    n.chr <- length(map1)
    maxloc <- max(c(unlist(lapply(map1,max)),unlist(lapply(map2,max))))

    if(!horizontal) {
      old.xpd <- par("xpd")
      old.xaxt <- par("xaxt")
      old.las <- par("las")
      par(xpd=TRUE,xaxt="n",las=1)
      on.exit(par(xpd=old.xpd,xaxt=old.xaxt,las=old.las))

      plot(0,0,type="n",ylim=c(maxloc,0),xlim=c(0.5,n.chr+0.5),
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
  invisible()
}


plot.cross <-
function (x, auto.layout = TRUE, pheno, ...) 
{
  old.yaxt <- par("yaxt")
  old.mfrow <- par("mfrow")
  on.exit(par(yaxt = old.yaxt, mfrow = old.mfrow))

  n.phe <- nphe(x)
  if(missing(pheno)) pheno <- 1:n.phe
  n.plot = length(pheno) + 2

  # automatically choose row/column structure for the plots
  if(auto.layout) {
    nr <- ceiling(sqrt(n.plot))
    nc <- ceiling((n.plot)/nr)
    par(mfrow = c(nr, nc))
  }

  plot.missing(x)
  plot.map(x)

#  if( is.numeric(pheno) )
#    pheno = names(x$pheno)[pheno]
  if(!is.numeric(pheno)) {
    temp <- match(pheno, names(x$pheno))
    if(any(is.na(temp))) 
      warning("Some phenotypes not found:",
              paste(pheno[is.na(temp)], collapse=" "))
    pheno <- temp[!is.na(temp)]
  }

  for(i in pheno) {
    if(!is.numeric(x$pheno[[i]])) {
      par(yaxt = "s")
      barplot(c(table(x$pheno[[i]])), axes = FALSE, xlab = paste("phe", i),
              ylab = "", main = colnames(x$pheno)[i], col = "white")
    }
    else hist(x$pheno[[i]], breaks = round(sqrt(nrow(x$pheno)) + 5),
              xlab = paste("phe", i), prob = TRUE, ylab = "", yaxt = "n",
              main = colnames(x$pheno)[i])
  }
  invisible()
}


##################################################r####################
#
# plot.geno: Plot genotypes for a specified chromosome, with likely
#           genotyping errors indicated. 
#
######################################################################

plot.geno <-
function(x, chr, ind, include.xo=TRUE, horizontal=TRUE,
         cutoff=3.5, min.sep=2, cex=1.2, ...)
{
  cross <- x  
  cross <- subset(cross,chr=chr)
  if(!missing(ind)) cross <- subset(cross, ind=ind)
  type <- class(cross)[1]
  
  if(type != "bc" && type != "f2" && type != "riself" && type != "risib")
    stop("Only available for backcross, intercross or RI strains.")

  if(is.na(match("errorlod",names(cross$geno[[1]])))) {
    warning("First running calc.errorlod.")
    cross <- calc.errorlod(cross,error.prob=0.01)
  }
  
  # indicators for apparent errors
  errors <- matrix(0,ncol=ncol(cross$geno[[1]]$data),
                   nrow=nrow(cross$geno[[1]]$data))
  dimnames(errors) <- dimnames(cross$geno[[1]]$data)

  top <- top.errorlod(cross,1,cutoff,FALSE)
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
  cross$geno[[1]]$map <- map
  
  data <- cross$geno[[1]]$data
  n.ind <- nrow(errors)

  color <- c("white","gray60","black","green","orange","red")

  if(include.xo) { # find crossover locations
    xoloc <- locate.xo(cross)
    xoloc <- data.frame(ind=rep(1:length(xoloc),sapply(xoloc,length)),
                        loc=unlist(xoloc))
  }

  if(horizontal==TRUE) {
    plot(0,0,type="n",xlab="Position (cM)",ylab="Individual",
         main=paste("Chromosome",names(cross$geno)[1]),
         ylim=c(0.5,n.ind+0.5),xlim=c(0,max(map)))
    segments(0,1:n.ind,max(map),1:n.ind)

    # AA genotypes
    tind <- rep(1:n.ind,length(map));tind[is.na(data)] <- NA
    ind <- tind; ind[!is.na(data) & data!=1] <- NA
    x <- rep(map,rep(n.ind,length(map)))
    points(x,ind,pch=16,col=color[1],cex=cex)
    points(x,ind,pch=1,cex=cex)

    # AB genotypes
    ind <- tind; ind[!is.na(data) & data!=2] <- NA
    if(type=="f2") {
      points(x,ind,pch=16,col=color[2],cex=cex)
      points(x,ind,pch=1,cex=cex)
    }
    else points(x,ind,pch=16,col=color[3],cex=cex) 

    if(type=="f2") {
      # BB genotypes
      ind <- tind; ind[!is.na(data) & data!=3] <- NA
      points(x,ind,pch=16,col=color[3],cex=cex)

      # not BB (D in mapmaker/qtl) genotypes
      ind <- tind; ind[!is.na(data) & data!=4] <- NA
      points(x,ind,pch=16,col=color[4],cex=cex)
      points(x,ind,pch=1,cex=cex)

      # not AA (C in mapmaker/qtl) genotypes
      ind <- tind; ind[!is.na(data) & data!=5] <- NA
      points(x,ind,pch=16,col=color[5],cex=cex)
      points(x,ind,pch=1,cex=cex)
    }

    # plot map
    u <- par("usr")
    segments(map,u[3],map,(u[3]+1)/2)
    segments(map,u[4],map,(n.ind+u[4])/2)

    if(any(errors)) {
      ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
      points(x,ind,pch=0,col=color[6],cex=cex+0.4,lwd=2)
    }

    if(include.xo) points(xoloc$loc,xoloc$ind,pch=4,col="blue",lwd=2)
  }
  else {
    plot(0,0,type="n",ylab="Position (cM)",xlab="Individual",
         main=paste("Chromosome",names(cross$geno)[1]),
         xlim=c(0.5,n.ind+0.5),ylim=c(max(map),0))
    segments(1:n.ind,0,1:n.ind,max(map))
    
    # AA genotypes
    tind <- rep(1:n.ind,length(map));tind[is.na(data)] <- NA
    ind <- tind; ind[!is.na(data) & data!=1] <- NA
    y <- rep(map,rep(n.ind,length(map)))
    points(ind,y,pch=16,col="white",cex=cex)
    points(ind,y,pch=1,cex=cex)

    # AB genotypes
    ind <- tind; ind[!is.na(data) & data!=2] <- NA
    if(type=="f2") {
      points(ind,y,pch=16,col=color[2],cex=cex)
      points(ind,y,pch=1,cex=cex)
    }
    else points(ind,y,pch=16,col=color[3],cex=cex)

    if(type=="f2") {
      # BB genotypes
      ind <- tind; ind[!is.na(data) & data!=3] <- NA
      points(ind,y,pch=16,col=color[3],cex=cex)

      # not BB genotypes
      ind <- tind; ind[!is.na(data) & data!=4] <- NA
      points(ind,y,pch=16,col=color[4],cex=cex)
      points(ind,y,pch=1,cex=cex)

      # not AA genotypes
      ind <- tind; ind[!is.na(data) & data!=5] <- NA
      points(ind,y,pch=16,col=color[5],cex=cex)
      points(ind,y,pch=1,cex=cex)
    }

    # plot map
    u <- par("usr")
    segments(u[1],map,(u[1]+1)/2,map)
    segments(u[2],map,(n.ind+u[2])/2,map)

    if(any(errors)) {
      ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
      points(ind,y,pch=0,col=color[6],cex=cex+0.4,lwd=2)
    }

    if(include.xo) points(xoloc$ind,xoloc$loc,pch=4,col="blue",lwd=2)

  }
  invisible()
}
    
######################################################################
#
# plot.info: Plot the proportion of missing information in the
#            genotype data.
#
######################################################################
plot.info <-
function(x,chr,method=c("both","entropy","variance"),...)
{
  cross <- x
  method <- match(match.arg(method),c("entropy","variance","both"))-1

  if(!missing(chr)) cross <- subset(cross,chr=chr)

  n.chr <- nchr(cross)
  results <- NULL

  if(is.na(match("prob",names(cross$geno[[1]])))) { # need to run calc.genoprob
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

#  gap <- attr(cross$geno[[1]]$prob,"off.end")*2+10 # gap between chr in plot
  gap <- 25

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
               as.integer(method),
               PACKAGE="qtl")

    if(method != 1) { # rescale entropy version
      if(n.gen==3) maxent <- 1.5*log(2)
      else maxent <- log(n.gen)
      info$info1 <- -info$info1/maxent
    }
    if(method != 0) { # rescale variance version
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
    if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
      w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
    rownames(z) <- w
    results <- rbind(results, z)
  }

  # check whether gap was included as an argument
  args <- list(...)
  if(is.na(match("gap",names(args)))) {
    if(method==0)
      plot.scanone(results,ylim=c(0,1),gap=gap,
                   main="Missing information",...)
    else if(method==1)
      plot.scanone(results,lodcolumn=2,ylim=c(0,1),gap=gap,
                   main="Missing information",...)
    else if(method==2)
      plot.scanone(results,results,lodcolumn=1:2,ylim=c(0,1),gap=gap,
                   main="Missing information",...)
  }
  else { # gap was included in ...
    if(method==0)
      plot.scanone(results,ylim=c(0,1),
                   main="Missing information",...)
    else if(method==1)
      plot.scanone(results,lodcolumn=2,ylim=c(0,1),
                   main="Missing information",...)
    else if(method==2)
      plot.scanone(results,results,lodcolumn=1:2,ylim=c(0,1),
                   main="Missing information",...)
  }

  colnames(results)[3:4] <- c("misinfo.entropy","misinfo.variance")

  class(results) <- c("scanone","data.frame")
  invisible(results)
}


# plot phenotypes against one or more markers
plot.pxg <-
function(x, marker, pheno.col = 1, jitter = 1, infer = TRUE, 
         pch, ylab, main, ...) 
{
  cross <- x
  type <- class(cross)[1]

  if(missing(pch)) pch <- par("pch")
  if(missing(ylab)) ylab <-  colnames(cross$pheno)[pheno.col] 

  oldlas <- par("las")
  on.exit(par(las = oldlas))
  par(las = 1)

  # find chromosomes containing the markers
  o <- sapply(cross$geno, function(a, b) !is.na(match(b, colnames(a$data))), 
              marker)
  if(length(marker)==1) o <- matrix(o,nrow=1)
  if(!all(apply(o,1,any))) {
    oo <- apply(o,1,any)
    err <- paste("Marker", marker[!oo], "not found")
    stop(err)
  }
  n.mark <- length(marker)
  o <- apply(o, 1, which)
  chr <- names(cross$geno)[o]
  uchr <- unique(chr)

  cross <- subset(cross, chr=uchr)
  map <- pull.map(cross)
  pos <- NULL
  for(i in seq(length(chr))) pos[i] <- map[[chr[i]]][marker[i]]
  chrtype <- sapply(cross$geno, class)
  if(length(chr) != length(uchr))
    chrtype <- chrtype[chr]
  
  # if X chromosome and backcross or intercross, get sex/direction data
  if(any(chrtype == "X") && (type == "bc" || type == "f2"))
    sexpgm <- getsex(cross)
  else sexpgm <- NULL

  # number of possible genotypes
  gen.names <- list()
  for(i in seq(length(chr))) 
    gen.names[[i]] <- getgenonames(type, chrtype[i], "full", sexpgm)
  n.gen <- sapply(gen.names, length)


  jitter <- jitter/10
  if(any(n.gen == 2)) jitter <- jitter * 0.75

  # function to determine whether genotype is fully known
  tempf <-
    function(x, type)
      {
        tmp <- is.na(x)
        if(type=="f2") tmp[!is.na(x) & x>3] <- TRUE
        if(type=="4way") tmp[!is.na(x) & x>4] <- TRUE
        tmp
      }

  # if infer=TRUE, fill in genotype data by a single imputation
  if(infer) {
    which.missing <- tempf(cross$geno[[chr[1]]]$data[, marker[1]],type)
    if(n.mark > 1) 
      for(i in 2:n.mark)
        which.missing <- which.missing | tempf(cross$geno[[chr[i]]]$data[,marker[i]],type)
    which.missing <- as.numeric(which.missing)

    cross <- fill.geno(cross, method = "imp")
  }
  else which.missing <- rep(1,nind(cross))

  # data to plot
  x <- cross$geno[[chr[1]]]$data[, marker[1]]
  if(n.mark > 1) 
    for(i in 2:n.mark)
      x <- cbind(x, cross$geno[[chr[i]]]$data[, marker[i]])
  else x <- as.matrix(x)
  y <- cross$pheno[, pheno.col]

  if(!infer) { # replace partially informative genotypes with NAs
    if(type == "f2") x[x > 3] <- NA
    if(type == "4way") x[x > 4] <- NA
  }

  # in case of X chromosome, recode some genotypes
  if(any(chrtype == "X") && (type == "bc" || type == "f2")) {
    ix = seq(n.mark)[chrtype == "X"]
    for(i in ix)
      x[, i] <- as.numeric(reviseXdata(type, "full", sexpgm,
                                       geno = as.matrix(x[, i])))
  }

  # save all of the data, returned invisibly
  data <- as.data.frame(x)
  names(data) <- marker
  for(i in marker) data[[i]] <- ordered(data[[i]])
  data$pheno <- y
  data$inferred <- which.missing

  # re-code the multi-marker genotypes
  if(n.mark > 1) {
    for(i in 2:n.mark)
      x[, 1] <- n.gen[i] * (x[, 1] - 1) + x[, i]
  }
  x <- x[, 1]

  # amount of jitter 
  u <- runif(nind(cross), -jitter, jitter)
  r <- (1 - 2 * jitter)/2

  # create plot
  plot(x + u, y, xlab = "Genotype", ylab = ylab, type = "n", 
       main = "", xlim = c(1 - r + jitter, prod(n.gen) + r + 
                    jitter), xaxt = "n")

  # marker names at top
  if(missing(main))
    mtext(paste(marker, collapse = "\n"), , 0.5, cex = max(2/n.mark, 
                                                   0.75))
  else
    title(main=main)
  
  abline(v = 1:prod(n.gen), col = "gray", lty = 3)

  if(length(pch) == 1) 
    pch = rep(pch, length(x))
  if(infer) {
    points((x + u)[which.missing == 1], y[which.missing == 
                     1], col = "red", pch = pch[which.missing == 1])
    points((x + u)[which.missing == 0], y[which.missing == 
                     0], pch = pch[which.missing == 0])
  }
  else points(x + u, y, pch = pch)
  sux = sort(unique(x))

  # add confidence intervals
  me <- se <- array(NA, prod(n.gen))
  me[sux] <- tapply(y, x, mean, na.rm = TRUE)
  se[sux] <- tapply(y, x, function(a) sd(a, na.rm = TRUE)/sqrt(sum(!is.na(a))))
  cols <- "blue"
  if(n.gen[n.mark] == 3) 
    cols <- c("blue", "purple", "red")
  else if(n.gen[n.mark] == 2) 
    cols <- c("blue", "red")
  segments(seq(prod(n.gen)) + jitter * 2, me, seq(prod(n.gen)) + 
           jitter * 4, me, lwd = 2, col = cols)
  segments(seq(prod(n.gen)) + jitter * 3, me - se, seq(prod(n.gen)) + 
           jitter * 3, me + se, lwd = 2, col = cols)
  segments(seq(prod(n.gen)) + jitter * 2.5, me - se, seq(prod(n.gen)) + 
           jitter * 3.5, me - se, lwd = 2, col = cols)
  segments(seq(prod(n.gen)) + jitter * 2.5, me + se, seq(prod(n.gen)) + 
           jitter * 3.5, me + se, lwd = 2, col = cols)

  # add genotypes below
  u <- par("usr")
  segments(1:prod(n.gen), u[3], 1:prod(n.gen), u[3] - diff(u[3:4]) * 
           0.015, xpd = TRUE)
  if(n.mark == 1) 
    tmp <- gen.names[[1]]
  else {
    tmp <- array(gen.names[[n.mark]], c(prod(n.gen), n.mark))
    for(i in (n.mark - 1):1) {
      tmpi <- rep(gen.names[[i]], rep(prod(n.gen[(i + 1):n.mark]), 
                                      n.gen[i]))
      if(i > 1) 
        tmpi <- rep(tmpi, prod(n.gen[1:(i - 1)]))
      tmp[, i] <- tmpi
    }
    tmp <- apply(tmp, 1, function(x) paste(x, collapse = "\n"))
  }
  text(1:prod(n.gen), u[3] - diff(u[3:4]) * 0.05, tmp, xpd = TRUE, 
       cex = max(0.5, 1.5/n.mark))

  invisible(data)

  # calculate return values?
#  if(any(which.missing == 0)) 
#    p.value <- anova(aov(y ~ x, subset = (which.missing == 
#                                          0)))[1, 5]
#  else p.value <- NA
#  names(p.value) <- NULL
#  tmp <- options(warn = -1)
#  form <- formula(paste("y ~", paste(marker, collapse = "*")))
#  if(any(is.na(me)) & n.mark > 2) {
#    formadd <- formula(paste("y ~", paste(marker, collapse = "+")))
#    fit <- aov(formadd, data, subset = (data$inferred == 
#                                        0))
#    full <- aov(form, data, subset = (data$inferred == 0))
#  }
#  else fit <- aov(form, data, subset = (data$inferred == 0))
#  tbl <- anova(fit, type = "marginal")
#  options(tmp)
#  p.value <- round(tbl$P[-nrow(tbl)], 4)
#  tmp = summary.lm(fit)
#  Rsq = tmp$r.sq
#  fstat = tmp$fstatistic
#  p.value = c(pf(fstat[1], fstat[2], fstat[3], lower = FALSE), 
#    p.value)
#  names(p.value) <- c("overall", dimnames(tbl)[[1]][-nrow(tbl)])
#  if(any(is.na(me)) & n.mark > 2) {
#    p.value["inter"] <- round(anova(fit, full)$P[2], 4)
#    fit = full
#  }
#  invisible(list(Rsq = Rsq, p.value = p.value, me = me, se = se, 
#                 fit = fit, data = data))

}

# end of plot.R
