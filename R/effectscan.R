######################################################################
#
# effectscan.R
#
# copyright (c) 2003-4, Hao Wu, The Jackson Laboratory
#                    with modifications by Karl W. Broman
# last modified Sep, 2004
# first written Jan, 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: effectscan
#
######################################################################

effectscan <-
function(cross, pheno.col=1, chr, ylim, gap=25,
         col=c("black","blue","red"), lty=c(1,2,3), lwd=2,
         mtick=c("line", "triangle"), main, add.legend=TRUE,
         ...)
{
  mtick <- match.arg(mtick)

  if(!missing(chr)) cross <- subset(cross,chr=chr)

  # remove individuals with missing phenotype (added by Karl)
  cross <- subset(cross, ind=!is.na(cross$pheno[,pheno.col]))

  pheno <- cross$pheno[,pheno.col]
  
  # loop thru all markers on all chromosomes
  # chromosome number
  chr <- NULL
  # x axis value for plot
  xvalue <- NULL
  # x-axis ticks and tick labels
  xtick <- NULL
  xticklabel <- NULL
  # y axis values, there're additive effect for bc
  # and additive and dominace effects for f2
  addeff <- NULL
  domeff <- NULL
  
  for(i in 1:nchr(cross)) {
    if(i != 1) tmp <- cross$geno[[i]]$map + max(xvalue) + gap
    else tmp <- cross$geno[[i]]$map
    xvalue <- c(xvalue, tmp)
    chr <- c(chr, rep(i, length(tmp)))
    xtick <- c(xtick, mean(c(min(tmp), max(tmp))))
    xticklabel <- c(xticklabel, names(cross$geno)[i])
    
    # find the y axis value in plot
    for(j in 1:dim(cross$geno[[i]]$data)[2]) {
      # the genotype for this marker
      geno <- cross$geno[[i]]$data[,j]
      if(class(cross)[1]=="bc") {
        # if this is back cross, 1 is A, 2 is H
        idx.1 <- which(geno==1)
        idx.2 <- which(geno==2)
        addeff <- c(addeff, mean(pheno[idx.1])-mean(pheno[idx.2]))
      }
      else if(class(cross)[1]=="f2"){
        # if this is F1, 1 is AA, 2 is AB, 3 is BB
        idx.1 <- which(geno==1)
        idx.2 <- which(geno==2)
        idx.3 <- which(geno==3)
        if(names(cross$geno[i]) != "X") {
          # automosomes 
          addeff <- c(addeff, mean(pheno[idx.1])-mean(pheno[idx.3]))
          # there's no dominance effect for X chromosome (only 2 genotypes)
          domeff <- c(domeff, mean(pheno[idx.2]) -
                      (mean(pheno[idx.1])+mean(pheno[idx.3]))/2)
        }
        else {
          # X chromosome, only 2 genotypes
          # there's no dominance effect, only additive effect
          addeff <- c(addeff, mean(pheno[idx.1])-mean(pheno[idx.2]))
        }
      }
      else {
        # other cross, implement later
      }
    }
  }

  # plot it
  # graphics parameters
  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=FALSE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))
  # line type, width, color
  if(length(lty)==1) lty <- rep(lty,3)
  if(length(lwd)==1) lwd <- rep(lwd,3)
  if(length(col)==1) col <- rep(col,3)

  if(missing(ylim)) {
    tmp <- c(addeff, domeff)
    tmp <- tmp[!is.na(tmp)]
    ylim <- range(tmp)
  }
  if(missing(main)) main <- "Effect scan plot"
  plot(0, 0, ylim=ylim, xlim=range(xvalue),type="n", xaxt="n", xlab="",
       ylab="", main=main, ...)
  for(i in 1:nchr(cross)) {
    # draw additive effects for this chromosome
    lines(xvalue[chr==i], addeff[chr==i], lwd=lwd[1],
          lty=lty[1], col=col[1])
    # if cross type is "f2" or others, add lines
    if(class(cross)[1]=="f2") {
      if(names(cross$geno[i]) != "X")
        lines(xvalue[chr==i], domeff[chr==i], lwd=lwd[2],
              lty=lty[2], col=col[2])
    }
  }

  # draw x axis ticks
  if(nchr(cross)>1) {
    axis(1, at=xtick, labels=xticklabel)
  }
  else {
    axis(1)
    title(xlab="Map position (cM)")
  }
  
  # add tick marker
  a <- par("usr")
  if(mtick=="line")
    rug(xvalue, 0.02, quiet=TRUE)
  else
    points(xvalue, rep(a[3]+diff(a[3:4])*0.04, length(xvalue)), pch=17, cex=1.5)

  # add legend (if requested and there are more than 2 lines)
  if(add.legend & !is.null(domeff)) {
    a <- par("usr")
    x.leg <- 0.15*a[1]+0.85*a[2]
    y.leg <- 0.05*a[3]+0.95*a[4]
    leg <- c("Additive", "Dominance")
    legend(x.leg, y.leg, leg, lty=lty[1:2], 
           col=col[1:2], cex=1, xjust=0.5)
  }

}

# end of effectscan.R
