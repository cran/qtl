#####################################################################
#
# arithscan.R
#
# copyright (c) 2005, Karl W Broman, Johns Hopkins University
# last modified Apr, 2005
# first written Mar, 2005
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Hao Wu (The Jackson Lab) wrote the imputation method
#
# Part of the R/qtl package
# Contains: +.scanone, -.scanone
#
######################################################################


"-.scanone" <-
function(e1,e2)
{
  if(missing(e2)) {
    class(e1) <- "data.frame"
    e1[,-(1:2)] <- -e1[,-(1:2)]
    class(e1) <- c("scanone","data.frame")
    return(e1)
  }

  class(e1) <- class(e2) <- "data.frame"
  if(nrow(e1) != nrow(e2)) {
    u1 <- levels(e1[,1])
    u2 <- levels(e2[,2])
    u <- unique(c(u1,u2))
    if(length(u) == 0) stop("Can't subtract; no chromosomes in common.")
    e1 <- e1[!is.na(match(e1[,1],u)),]
    e2 <- e2[!is.na(match(e2[,1],u)),]
    if(nrow(e1) != nrow(e2) || any(e1[,1] != e2[,1]) || max(abs(e1[,2]-e2[,2])) < 0.01)
      stop("Can't subtract; arguments not compatible")
  }
  nc1 <- ncol(e1)
  nc2 <- ncol(e2)
  nc <- min(c(nc1,nc2))
  e1 <- e1[,1:nc]
  e1[,3:nc] <- e1[,3:nc] - e2[,3:nc]
  class(e1) <- c("scanone","data.frame")
  e1
}

"+.scanone" <-
function(e1,e2)
{
  if(missing(e2)) return(e1)

  class(e1) <- class(e2) <- "data.frame"
  if(nrow(e1) != nrow(e2)) {
    u1 <- levels(e1[,1])
    u2 <- levels(e2[,2])
    u <- unique(c(u1,u2))
    if(length(u) == 0) stop("Can't subtract; no chromosomes in common.")
    e1 <- e1[!is.na(match(e1[,1],u)),]
    e2 <- e2[!is.na(match(e2[,1],u)),]
    if(nrow(e1) != nrow(e2) || any(e1[,1] != e2[,1]) || max(abs(e1[,2]-e2[,2])) < 0.01)
      stop("Can't subtract; arguments not compatible")
  }
  nc1 <- ncol(e1)
  nc2 <- ncol(e2)
  nc <- min(c(nc1,nc2))
  e1 <- e1[,1:nc]
  e1[,3:nc] <- e1[,3:nc] + e2[,3:nc]
  class(e1) <- c("scanone","data.frame")
  e1
}

# end of arithscan.R
