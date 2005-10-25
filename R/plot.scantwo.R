######################################################################
#
# plot.scantwo.R
#
# copyright (c) 2001-5, Karl W Broman, Johns Hopkins University,
#                       Hao Wu and Brian Yandell
# last modified Oct, 2005
# first written Nov, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Hao Wu (The Jackson Lab) wrote the initial code
#
# Part of the R/qtl package
# Contains: plot.scantwo, subset.scantwo
#
######################################################################

plot.scantwo <-
function(x, chr, incl.markers = FALSE, zlim,
         lower = c("joint", "add", "cond-int", "cond-add"), nodiag = TRUE,
         contours = FALSE, main, zscale = TRUE, point.at.max=FALSE,
         col.scheme = c("redblue","cm","gray","heat","terrain","topo"),
         gamma = 1, ...)
{
  col.scheme <- match.arg(col.scheme)

  if(!missing(chr)) 
    x <- subset(x, chr=chr)
  chr <- as.character(unique(x$map[,1]))

  lower <- match.arg(lower)
  if(!any(class(x) == "scantwo")) 
    stop("Input variable is not an object of class scantwo!")
  lod <- x$lod
  map <- x$map

  # backward compatibility for previous version of R/qtl
  if(is.na(match("scanoneX",names(x)))) {
    warning("It would be best to re-run scantwo() with the R/qtl version 0.98 or later.")
    scanoneX <- NULL
  }
  else scanoneX <- x$scanoneX

  # if incl.markers is FALSE, drop positions
  #     for which third column of map is 0
  if(!incl.markers && any(map[, 3] == 0)) {
    o <- (map[, 3] == 1)
    lod <- lod[o, o]
    map <- map[o, ]
    if(!is.null(scanoneX)) scanoneX <- scanoneX[o]
  }

  if(all(diag(lod) < 1e-14) && lower != "joint") 
    stop("Need to run scantwo with run.scanone=TRUE.")

  # pull out single-QTL LODs
  if(lower=="cond-int" || lower=="cond-add") {
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
    q1[q2>q1] <- q2[q2>q1]
  }

  # replace joint LOD with LOD[q1,q2] - max{LOD[q1],LOD[q2]}
  if(lower == "cond-int") {
    lod[lower.tri(lod)] <- (lod - q1)[lower.tri(lod)] 
  }
  else if(lower == "cond-add") {
    lod[lower.tri(lod)] <- (lod - t(lod) - q1)[lower.tri(lod)]
  }
  else if(lower == "add") {
    lod[lower.tri(lod)] <- (lod - t(lod))[lower.tri(lod)]
  }
  
  if(nodiag) diag(lod) <- 0

  # deal with bad LOD score values
  if(any(is.na(lod) | lod < -1e-06 | lod == Inf)) {
    warning("Some LOD scores NA, Inf or < 0; set to 0")
    lod[is.na(lod) | lod < 0 | lod == Inf] <- 0
  }

  if(missing(zlim)) { # no given zlim
    # calculate the zlim for interactive and joint
    zlim.int <- max(lod[row(lod) < col(lod)])
    zlim.jnt <- max(lod[row(lod) >= col(lod)])
  }
  else {
    zlim.int <- zlim[2]
    zlim.jnt <- zlim[1]
  }

  # rescale the data in upper triangle based on zlims.jnt
  lod[row(lod) < col(lod)] <- lod[row(lod) < col(lod)] * zlim.jnt/zlim.int
  if(missing(zlim)) 
    zlim.jnt <- max(lod)

  # make sure LOD values are below (0,zlim.jnt) or update zlim.jnt
  if(max(lod) > zlim.jnt) {
    warning("LOD values out of range; updating zlim.")
    temp <- max(lod)
    zlim.int <- zlim.int * temp/zlim.jnt
    zlim.jnt <- temp
  }

  # save old par parameters, to restore them on exit
  old.mar <- par("mar")
  old.las <- par("las")
  old.mfrow <- par("mfrow")
  on.exit(par(las = old.las, mar = old.mar, mfrow = old.mfrow))
  par(las = 1)
  if(zscale) {
    layout(cbind(1, 2), c(6, 1))
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
  if( gamma < 0 && col.scheme == "redblue")
    stop( "gamma must be non-negative" )
  cols <- switch(col.scheme,
                 gray = if( gamma <= 0) rev(gray(seq(0,1,len=256)))
                   else rev(gray(log(seq(1,exp(gamma),len=256))/gamma)),
                 heat = heat.colors(256),
                 terrain = terrain.colors(256),
                 topo = topo.colors(256),
                 cm = cm.colors(256),
                 redblue = rev(rainbow(256, start = 0, end = 2/3,gamma=gamma)))

  if(length(chr) > 1)
    image(1:ncol(lod), 1:nrow(lod), lod, ylab = "Chromosome", 
          xlab = "Chromosome", zlim = c(0, zlim.jnt), col = cols,
          xaxt = "n", yaxt = "n")
  else
    image(map[,2], map[,2], lod, ylab = "Position (cM)", 
          xlab = "Position (cM)", zlim = c(0, zlim.jnt), col = cols)
          

  # plot point at maximum, if requested
  if(point.at.max) {
    temp <- lod
    temp[upper.tri(temp)] <- -50
    temp[diag(temp)] <- -50
    wh <- which(temp == max(temp), arr.ind=TRUE)
    if(length(chr) > 1)
      points(wh,rev(wh),pch=4,lwd=2)
    else {
      points(map[wh[,1],2],map[wh[,2],2],pch=4,lwd=2,col="blue")
      points(map[wh[,2],2],map[wh[,1],2],pch=4,lwd=2,col="blue")
    }
  }

  # add contours if requested
  if(any(contours) > 0) {
    if(is.logical(contours))
      contours = 1.5
    tmp = lod
    tmp[row(lod) < col(lod)] <- NA
    if(length(chr) > 1) 
      thepos <- 1:ncol(lod)
    else thepos <- map[,2]
    contour(thepos, thepos, tmp, add = TRUE,drawlabels=FALSE,
            levels = max(tmp,na.rm=TRUE) - contours, col = "blue",
            lwd = 2)
    tmp = lod
    tmp[row(lod) > col(lod)] <- NA
    contour(thepos, thepos, tmp, add = TRUE,drawlabels=FALSE,
            levels = max(tmp,na.rm=TRUE) - contours * zlim.jnt/zlim.int,
            col = "blue", lwd = 2)
  }

  if(length(chr) > 1) {
    # calculate how many markers in each chromesome
    n.mar <- NULL
    for(i in 1:length(chr)) n.mar[i] <- sum(map[, 1] == chr[i])

    # plot lines at the chromosome boundaries
    if(length(chr) > 1)
      wh <- c(0.5, cumsum(n.mar) + 0.5)
    abline(v = wh, xpd = FALSE)
    abline(h = wh, xpd = FALSE)

    # add chromesome numbers
    a <- par("usr")
    for(i in 1:length(n.mar)) {
      text(mean(wh[i + c(0, 1)]), a[3] - diff(a[3:4]) * 0.025, 
           chr[i], xpd = TRUE, adj = c(0.5, 1))
      segments(mean(wh[i + c(0, 1)]), a[3],
               mean(wh[i + c(0, 1)]), a[3] - diff(a[3:4]) * 0.01, xpd = TRUE)
      text(a[1] - diff(a[1:2]) * 0.025, mean(wh[i + c(0, 1)]), 
           chr[i], xpd = TRUE, adj = c(1, 0.5))
      segments(a[1], mean(wh[i + c(0, 1)]), a[1] - diff(a[1:2]) * 
               0.01, mean(wh[i + c(0, 1)]), xpd = TRUE)
    }
  }

  # add title
  if(!missing(main)) 
    title(main = main)

  if(zscale) {
    # plot the colormap
    par(mar = c(5, 2, 4, 2) + 0.1)
    colorstep <- zlim.jnt/255
    image(x = 1:1, y = seq(0, zlim.jnt, colorstep), z = matrix(c(1:256), 1, 256),
          zlim = c(1, 256), ylab = "", xlab = "", 
          xaxt = "n", yaxt = "n", col = cols)

    # make sure there's a box around it
    u <- par("usr")
    abline(v = u[1:2], xpd = FALSE)
    abline(h = u[3:4], xpd = FALSE)
    if(any(contours) > 0) {
      for(i in seq(length(contours))) {
        segments(mean(u[1:2]),
                 max(lod[row(lod) > col(lod)]) - contours[i],
                 u[2], max(lod[row(lod) > col(lod)]) - contours[i], 
                 xpd = FALSE, col = "blue", lwd = 2)
        segments(u[1], max(lod[row(lod) < col(lod)]) - contours[i] * 
                 zlim.jnt/zlim.int, mean(u[1:2]),
                 max(lod[row(lod) < col(lod)]) - contours[i] * zlim.jnt / zlim.int,
                 xpd = FALSE, col = "blue", lwd = 2)
      }
    }

    # figure out how big the axis labels should be
    fin <- par("fin")[1] # figure width in inches
    pin <- par("pin")[1] # plot width in inches
    mai <- par("mai")[2] # margin width in inches
                         # note: pin + 2*mai = fin
    xlen.mar <- mai/pin * diff(u[1:2])

    # axis for joint LODs
    yloc <- pretty(c(0, zlim.jnt), 4)
    yloc <- yloc[yloc <= u[4]]
    segments(u[2], yloc, u[2] + xlen.mar/4, yloc, xpd = TRUE)
    text(u[2] + xlen.mar/3, yloc, as.character(yloc), xpd = TRUE, adj = 0)

    # axis for int've LODs
    yloc <- pretty(c(0, zlim.int), 4)
    yloc.rev <- yloc * zlim.jnt/zlim.int
    yloc <- yloc[yloc.rev <= u[4]]
    yloc.rev <- yloc.rev[yloc.rev <= u[4]]
    segments(u[1], yloc.rev, u[1] - xlen.mar/4, yloc.rev, xpd = TRUE)
    text(u[1] - xlen.mar/3, yloc.rev, as.character(yloc), xpd = TRUE, adj = 1)
  }

  invisible()
}

######################################################################
#
# subset.scantwo
#
######################################################################

#subset.scantwo <-
#function(x, chr, ...)   
#{
#  if(missing(chr) || length(chr) == 0) return(x)
#
#  a <- unique(x$map[,1])
#  if(is.numeric(chr) && all(chr < 0)) 
#    chr <- a[chr]
#  else chr <- a[match(chr,a)]
#
#  newgroups <- groups <- vector("list",length(chr))
#  curmax <- 0
#  for(i in 1:length(chr)) {
#    groups[[i]] <- which(x$map[,1]==chr[i])
#    newgroups[[i]] <- 1:length(groups[[i]]) + curmax
#    curmax <- curmax + length(groups[[i]])
#  }
#
#  g <- unlist(groups)
#  x$map <- x$map[g,]
#
#  lod <- matrix(ncol=length(g),nrow=length(g))
#  for(i in 1:length(chr)) {
#    lod[newgroups[[i]],newgroups[[i]]] <- x$lod[groups[[i]],groups[[i]]]
#    if(i < length(chr))
#      for(j in (i+1):length(chr)) {
#        if(groups[[i]][1] < groups[[j]][2]) {
#          lod[newgroups[[i]],newgroups[[j]]] <- x$lod[groups[[i]],groups[[j]]]
#          lod[newgroups[[j]],newgroups[[i]]] <- x$lod[groups[[j]],groups[[i]]]
#        }
#        else {
#          lod[newgroups[[j]],newgroups[[i]]] <- t(x$lod[groups[[i]],groups[[j]]])
#          lod[newgroups[[i]],newgroups[[j]]] <- t(x$lod[groups[[j]],groups[[i]]])
#        }
#      }
#  }
#  x$lod <- lod
#  x
#}

subset.scantwo <-
function(x, chr, ...)
{
  if(missing(chr) || length(chr)==0) return(x)

  a <- unique(x$map[,1])
  if(is.numeric(chr) && all(chr < 0)) 
    chr <- a[chr]
  else chr <- a[match(chr,a)]

  wh <- !is.na(match(x$map[,1],chr))

  if(length(wh) == 0) return(x)

  x$map <- x$map[wh,]
  x$lod <- x$lod[wh,wh]
  if(!is.null(x$scanoneX))
    x$scanoneX <- x$scanoneX[wh]

  x
}


# end of plot.scantwo.R
