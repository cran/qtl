######################################################################
#
# effectplot.R
#
# copyright (c) 2002-4, Hao Wu, The Jackson Laboratory
#                     and Karl W. Broman, Johns Hopkins University
# last modified Jul, 2004
# first written Jul, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: effectplot
#
######################################################################

effectplot <-
function (cross, pheno.col = 1, mname1, mark1, geno1, mname2, 
          mark2, geno2, main, ylim, add.legend = TRUE) 
{
  if(!sum(class(cross) == "cross")) 
    stop("The first input variable must be  an object of class cross")
  if(pheno.col > nphe(cross)) 
    stop("Input pheno.col is wrong")

  # local variables
  n.ind <- nind(cross)
  pheno <- cross$pheno[, pheno.col]
  type <- class(cross)[1]
  chrtype1 <- chrtype2 <- "A"
  gennames1 <- gennames2 <- NULL

  # Get marker 1 genotype data
  if(missing(mark1)) { # no data given
    if(missing(mname1)) 
      stop("Either mname1 or mark1 must be specified.")

    # find chromosome containing marker
    o <- sapply(cross$geno, function(a, b)
                !is.na(match(b, colnames(a$data))), mname1)

    if(!any(o)) {
      err <- paste("Marker", mname1, "not found")
      stop(err)
    }
    chr1 <- names(cross$geno)[o]
    chrtype1 <- class(cross$geno[[chr1]])

    # get genotype data
    if(!any(colnames(cross$geno[[chr1]]$data) == mname1)) {
      err <- paste("Marker", mname1, "not found.")
      stop(err)
    }
    mark1 <- cross$geno[[chr1]]$data[, mname1]

    # if X chr and backcross or intercross, get sex/dir data + revise data
    if(chrtype1 == "X" && (type == "bc" || type == "f2" || type == "f2ss")) {
      sexpgm <- getsex(cross)
      mark1 <- as.numeric(reviseXdata(type, "standard", sexpgm, 
                                      geno = as.matrix(mark1)))
      gennames1 <- getgenonames(type, chrtype1, "standard", sexpgm)
    }
  }
  else {
    if(length(mark1) != n.ind) 
      stop("Marker 1 data is the wrong length")
    if(missing(mname1)) 
      mname1 <- "Marker 1"
  }

  # Deal with marker 2
  if(!missing(mname2) || !missing(mark2)) {
    if(missing(mark2)) {

      # find chromosome containing marker
      o <- sapply(cross$geno, function(a, b)
                  !is.na(match(b, colnames(a$data))), mname2)
      if(!any(o)) {
        err <- paste("Marker", mname2, "not found")
        stop(err)
      }
      chr2 <- names(cross$geno)[o]
      chrtype2 <- class(cross$geno[[chr2]])

      # get genotype data
      if(!any(colnames(cross$geno[[chr2]]$data) == mname2)) {
        err <- paste("Marker", mname2, "not found.")
        stop(err)
      }
      mark2 <- cross$geno[[chr2]]$data[, mname2]

      # if X chr and backcross or intercross, get sex/dir data + revise data
      if(chrtype2 == "X" && (type == "bc" || type == "f2" || type == "f2ss")) {
        sexpgm <- getsex(cross)
        mark2 <- as.numeric(reviseXdata(type, "standard", sexpgm, 
                                        geno = as.matrix(mark2)))
        gennames2 <- getgenonames(type, chrtype2, "standard", sexpgm)
      }
    }
    else {
      if(length(mark2) != n.ind) 
        stop("Marker 2 data is the wrong length")
      if(missing(mname2)) 
        mname2 <- "Marker 2"
    }
  }
  else mark2 <- NULL

  # drop data for individuals with missing phenotypes or genotypes
  if(is.null(mark2)) {
    keepind <- !is.na(pheno) & !is.na(mark1)
    mark1 <- mark1[keepind]
    pheno <- pheno[keepind]
  }
  else {
    keepind <- !is.na(pheno) & !is.na(mark1) & !is.na(mark2)
    mark1 <- mark1[keepind]
    mark2 <- mark2[keepind]
    pheno <- pheno[keepind]
  }

  # adjust marker data and get level names
  if(!missing(geno1)) {
    if(length(unique(mark1)) > length(geno1)) 
      stop("geno1 is too short.")
    mark1 <- as.numeric(factor(mark1), levels = sort(unique(mark1)))
  }
  else {
    if(!is.null(gennames1)) 
      geno1 <- gennames1
    else if(is.factor(mark1)) {
      geno1 <- levels(mark1)
      mark1 <- as.numeric(mark1)
    }
    else {
      if(type == "bc") 
        geno1 <- c("AA", "AB")
      else if(type == "f2" || type == "f2ss") 
        geno1 <- c("AA", "AB", "BB")
      else if(type == "riself" || type == "risib") 
        geno1 <- c("AA", "BB")
      else if(type == "4way") 
        geno1 <- c("AC", "BC", "AD", "BD")
      if(length(unique(mark1)) > length(geno1)) 
        geno1 <- c(geno1, rep("?", length(unique(mark1)) - 
                              length(geno1)))
    }
  }
  if(!is.null(mark2)) {
    if(!missing(geno2)) {
      if(length(unique(mark2)) > length(geno2)) 
        stop("geno2 is too short.")
      mark2 <- as.numeric(factor(mark2), levels = sort(unique(mark2)))
    }
    else {
      if(!is.null(gennames2)) 
        geno2 <- gennames2
      else if(is.factor(mark2)) {
        geno2 <- levels(mark2)
        mark2 <- as.numeric(mark2)
      }
      else {
        if(type == "bc") 
          geno2 <- c("AA", "AB")
        else if(type == "f2" || type == "f2ss") 
          geno2 <- c("AA", "AB", "BB")
        else if(type == "riself" || type == "risib") 
          geno2 <- c("AA", "BB")
        else if(type == "4way") 
          geno2 <- c("AC", "BC", "AD", "BD")
        if(length(unique(mark2)) > length(geno2)) 
          geno2 <- c(geno2, rep("?", length(unique(mark2)) - 
                                length(geno2)))
      }
    }
  }
  ngen1 <- length(geno1)
  if(!is.null(mark2)) 
    ngen2 <- length(geno2)

  # calculate means and stds for interaction
  # and make output object
  # the output will be a data frame. For two-marker case,
  # the rows corresponding to the first marker and the columns
  # corresponding to the second marker
  result <- NULL
  if(is.null(mark2)) {
    means <- tapply(pheno, mark1, mean, na.rm = TRUE)
    ses <- tapply(pheno, mark1, function(a) sd(a, na.rm = TRUE)/sqrt(sum(!is.na(a))))
    lo <- means - ses
    hi <- means + ses
    # Note: rows are marker 1 and columns are marker 2

    if(length(means) != length(geno1)) {
      warning("Number of genotypes is different than length(geno1).")
      if(length(means) < length(geno1)) 
        geno1 <- geno1[1:length(means)]
      else geno1 <- c(geno1, rep("?", length(means) - length(geno1)))
      ngen1 <- length(geno1)
    }
    result$Means <- means
    names(result$Means) <- paste(mname1, geno1, sep = ".")
    result$SDs <- ses
    names(result$SDs) <- paste(mname1, geno1, sep = ".")
  }
  else {
    means <- tapply(pheno, list(mark1, mark2), mean, na.rm = TRUE)
    ses <- tapply(pheno, list(mark1, mark2), function(a) sd(a, 
                                                            na.rm = TRUE)/sqrt(sum(!is.na(a))))
    lo <- means - ses
    hi <- means + ses
    if(nrow(means) != length(geno1)) {
      warning("Number of genotypes in marker 1 is different than length(geno1).")
      if(nrow(means) < length(geno1)) 
        geno1 <- geno1[1:nrow(means)]
      else geno1 <- c(geno1, rep("?", nrow(means) - length(geno1)))
      ngen1 <- length(geno1)
    }
    if(ncol(means) != length(geno2)) {
      warning("Number of genotypes in marker 2 is different than length(geno2).")
      if(ncol(means) < length(geno2)) 
        geno2 <- geno2[1:ncol(means)]
      else geno2 <- c(geno2, rep("?", ncol(means) - length(geno2)))
      ngen2 <- length(geno2)
    }
    result$Means <- as.data.frame(means)
    rownames(result$Means) <- paste(mname1, geno1, sep = ".")
    colnames(result$Means) <- paste(mname2, geno2, sep = ".")
    result$SDs <- as.data.frame(ses)
    rownames(result$SDs) <- paste(mname1, geno1, sep = ".")
    colnames(result$SDs) <- paste(mname2, geno2, sep = ".")
  }

  ######### Draw the figure ############
  # graphics parameters
  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd = FALSE, las = 1)
  on.exit(par(xpd = old.xpd, las = old.las))

  # colors (for case of two markers)
  if(ngen1 <= 5) 
    int.color <- c("black", "red", "blue", "orange", "green")[1:ngen1]
  else int.color <- c("black", rainbow(ngen1 - 1, start = 0, 
                                       end = 2/3))
  # plot title
  if(missing(main)) {
    if(is.null(mark2)) 
      main <- paste("Effect plot for", mname1)
    else main <- paste("Interaction plot for", mname1, "and", 
                       mname2)
  }

  # y axis limits
  if(missing(ylim)) {
    ylimits <- range(c(lo, means, hi), na.rm = TRUE)
    ylimits[2] <- ylimits[2] + diff(ylimits) * 0.1
  }
  else ylimits <- ylim

  # x axis limits
  if(is.null(mark2)) { # one marker
    u <- sort(unique(mark1))
    d <- diff(u[1:2])
    xlimits <- c(min(mark1) - d/4, max(mark1) + d/4)
  }
  else { # two markers
    u <- sort(unique(mark2))
    d <- diff(u[1:2])
    xlimits <- c(min(mark2) - d/4, max(mark2) + d/4)
  }

  ## fix of x limits
  d <- 1
  xlimits <- c(1 - d/4, length(u) + d/4)

  if(is.null(mark2)) { # single marker
    # plot the means
    plot(1:ngen1, means, main = main, xlab = mname1, ylab = names(cross$pheno)[pheno.col], 
         pch = 1, col = "black", ylim = ylimits, xaxt = "n", 
         type = "b", xlim = xlimits)
    # confidence limits
    for(i in 1:ngen1) {
      if(!is.na(lo[i]) && !is.na(hi[i])) 
        lines(c(i, i), c(lo[i], hi[i]), pch = 3, col = "black", 
              type = "b", lty = 3)
    }

    # X-axis ticks
    a <- par("usr")
    ystart <- a[3]
    yend <- ystart - diff(a[3:4]) * 0.02
    ytext <- ystart - diff(a[3:4]) * 0.05
    for(i in 1:ngen1) {
      lines(x = c(i, i), y = c(ystart, yend), xpd = TRUE)
      text(i, ytext, geno1[i], xpd = TRUE)
    }
  }
  else { # two markers
    # plot the first genotype of marker 1
    plot(1:ngen2, means[1, ], main = main, xlab = mname2, 
         ylab = names(cross$pheno)[pheno.col], pch = 1, col = int.color[1], 
         ylim = ylimits, xaxt = "n", type = "b", xlim = xlimits)
    # confidence limits
    for(i in 1:ngen2) {
      if(!is.na(lo[1, i]) && !is.na(hi[1, i])) 
        lines(c(i, i), c(lo[1, i], hi[1, i]), pch = 3, 
              col = int.color[1], type = "b", lty = 3)
    }
    for(j in 2:ngen1) { # for the rest of genotypes for Marker 1
      lines(1:ngen2, means[j, ], col = int.color[j], pch = 1, 
            type = "b")
      # confidence limits
      for(i in 1:ngen2) {
        if(!is.na(lo[j, i]) && !is.na(hi[j, i])) 
          lines(c(i, i), c(lo[j, i], hi[j, i]), pch = 3, 
                col = int.color[j], type = "b", lty = 3)
      }
    }

    # draw X-axis ticks
    a <- par("usr")
    ystart <- a[3]
    yend <- ystart - diff(a[3:4]) * 0.02
    ytext <- ystart - diff(a[3:4]) * 0.05
    for(i in 1:ngen2) {
      lines(x = c(i, i), y = c(ystart, yend), xpd = TRUE)
      text(i, ytext, geno2[i], xpd = TRUE)
    }

    # add legend
    if(add.legend) {
      col <- int.color[1:ngen1]
      u <- sort(unique(mark2))
      x.leg <- mean(u[ngen2 - (0:1)])
      y.leg <- a[4] - diff(a[3:4]) * 0.05
      y.leg2 <- a[4] - diff(a[3:4]) * 0.03
      legend(x.leg, y.leg, geno1, lty = 1, pch = 1, col = col, 
             cex = 1, xjust = 0.5)
      text(x.leg, y.leg2, mname1)
    }
  }

  return(invisible(result))
}

# end of effectplot.R

