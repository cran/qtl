\name{plot.scantwo}
\alias{plot.scantwo}

\title{Plot LOD scores for a two-dimensional genome scan}

\description{
  Plot the joint and epistasis LOD scores for a two-dimensional genome
  scan with a two-QTL model. 
}

\usage{
plot.scantwo(x, chr, incl.markers=FALSE, zlim,
             contours=FALSE, main, zscale=TRUE, \dots)
}			    

\arguments{
  \item{x}{An object of class \code{"scantwo"}, as output by
    \code{link[qtl]{scantwo}}.} 
  \item{chr}{Vector specifying which chromosomes to plot.}
  \item{incl.markers}{If FALSE, plot LOD scores on an evenly
    spaced grid (not including the results at the markers).}
  \item{zlim}{A vector of length 2 (optional), indicating the z limits
    for the joint LOD and the epistasis LOD scores, respectively.  If
    missing, the maximum limits are used for each.}
  \item{contours}{If TRUE, add contours to the plot, using a call to
    \code{\link[base]{contour}}.}
  \item{main}{An optional title for the plot.}
  \item{zscale}{If TRUE, a color scale is plotted at the right.}
  \item{\dots}{Ignored at this point.}
}

\value{None.}

\details{
  Uses \code{\link[base]{image}} to plot a grid of LOD scores.  The
  upper triangle contains the epistasis LOD scores; the lower triangle
  contains the joint LOD scores.  The diagonal contains either all zeros
  or the main effects LOD scores.

  If a color scale is plotted (\code{zscale=TRUE}), the axis on the
  left indicates the scale for the epistasis LOD scores (the upper
  triangle), while the axis on the right indicates the joint LOD
  scores (the lower triangle).  Note that the axis labels get screwed up
  if you change the size of the figure window; you'll need to redo the
  plot. 

  Note that the numbers plotted with the contours (if
  \code{contours=TRUE}) are for the joint LOD scores, and may be
  incorrect for the epistasis LOD scores.  It appears to be difficult to
  fix this.   
}


\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2, step=10)
out.2dim <- scantwo(fake.f2, method="hk", chr=c(1,6,13))
plot(out.2dim)
plot(out.2dim, zlim=c(7.5,7.5))
}

\seealso{ \code{\link[qtl]{scantwo}},
  \code{\link[qtl]{summary.scantwo}}, \code{\link[qtl]{plot.scanone}} }

\author{Hao Wu, \email{hao@jax.org}; Karl W Broman,
  \email{kbroman@jhsph.edu} } 

\keyword{hplot}
