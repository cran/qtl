######################################################################
#
# scantwo.R
#
# copyright (c) 2001-2, Karl W Broman, Johns Hopkins University,
#                       and Hao Wu, The Jackson Lab.
# last modified Oct, 2002
# first written Nov, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Hao Wu (The Jackson Lab) wrote the initial code
#
# Part of the R/qtl package
# Contains: plot.scantwo
#
######################################################################

plot.scantwo <- 
function(x, chr, incl.markers=FALSE, zlim,
         lower=c("cond-int","cond-add","joint"), nodiag=TRUE,
         contours=FALSE, main, zscale=TRUE,...)
{
  lower <- match.arg(lower)
  if( !any(class(x) == "scantwo") )
    stop("Input variable is not an object of class scantwo!")
  
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

  if(all(diag(lod) < 1e-14))
    stop("Need to run scantwo with run.scanone=TRUE.")

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

  # replace joint LOD with LOD[q1,q2] - max{LOD[q1],LOD[q2]}
  if(lower=="cond-int") {
    for(i in 2:nrow(lod))
      for(j in 1:(i-1)) 
        lod[i,j] <- max(c(0,lod[i,j]-max(c(lod[i,i],lod[j,j]))))
  }
  else if(lower=="cond-add") {
    for(i in 2:nrow(lod))
      for(j in 1:(i-1))
        lod[i,j] <- max(c(0,lod[i,j]-lod[j,i]-max(c(lod[i,i],lod[j,j]))))
  }
    
  if(nodiag) diag(lod) <- 0

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

  image( 1:ncol(lod), 1:nrow(lod), lod, ylab="Chromosome", xlab="Chromosome",
         zlim=c(0,zlim.jnt), col=rev(rainbow(256,start=0,end=2/3)),
        xaxt="n", yaxt="n")

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
  for(i in 1:length(n.mar)) {
    text(mean(wh[i+c(0,1)]),a[3]-diff(a[3:4])*0.025,chr[i],xpd=TRUE,
         adj=c(0.5,1))
    segments(mean(wh[i+c(0,1)]), a[3], mean(wh[i+c(0,1)]),
             a[3]-diff(a[3:4])*0.01,xpd=TRUE)
    text(a[1]-diff(a[1:2])*0.025,mean(wh[i+c(0,1)]),chr[i],xpd=TRUE,
         adj=c(1,0.5))
    segments(a[1], mean(wh[i+c(0,1)]),
             a[1]-diff(a[1:2])*0.01, mean(wh[i+c(0,1)]), xpd=TRUE)

  }
  
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

