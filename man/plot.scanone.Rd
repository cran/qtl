\name{plot.scanone}
\alias{plot.scanone}

\title{Plot LOD curves}

\description{
  Plot the LOD curve for a genome scan with a single-QTL model (the
  output of \code{\link[qtl]{scanone}}).
}

\usage{
plot.scanone(x, x2, x3, chr, lodcolumn=3, incl.markers=TRUE, ylim,
             lty=c(1,2,3), col="black", lwd=2, add=FALSE,
             gap=25, main, \dots)
}			    

\arguments{
  \item{x}{A data.frame with at least three columns.  The first
    three columns should be the following (any further columns are
    ignored):} 
  \tabular{ll}{
     \code{chr} \tab Chromosome ID \cr
     \code{pos} \tab Position on genetic map \cr
     \code{lod}  \tab LOD score (or other statistic)}
  \item{x2}{data.frame for a second genome scan (optional)}
  \item{x3}{data.frame for a third genome scan (optional)}
  \item{chr}{Vector specifying which chromosomes to plot. (The
    chromosomes must be specified by name.)}
  \item{lodcolumn}{An integer, or vector of 3 integers, indicating which
    columns of the input should be plotted (generally column 3).}
  \item{incl.markers}{Indicate whether to plot line segments at the
  marker locations.}
  \item{ylim}{Limits for y-axis (optional).}
  \item{lty}{Line types; a vector of length 1 or 3.}
  \item{col}{Line colors; a vector of length 1 or 3.}
  \item{lwd}{Line widths; a vector of length 1 or 3.}
  \item{add}{If TRUE, add to a current plot.}
  \item{gap}{Gap separating chromosomes (in cM).}
  \item{main}{An optional title for the plot.}
  \item{\dots}{Passed to the function \code{\link[base]{plot}} when it
    is called.}
}
\value{None.}

\details{
  This function allows you to plot the results of up to three genome
  scans against one another.  Such objects must conform with each
  other. 

  One may alternatively use the argument \code{add} to add the plot of
  an additional genome scan to the current figure, but some care is
  required.  
}

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2,step=2.5)
out.mr <- scanone(fake.f2, method="mr")
out.em <- scanone(fake.f2, chr=c(1,13), method="em")
plot(out.mr)
plot(out.mr, out.em, chr=c(1,13))
plot(out.mr, out.em, , c(1,13), lty=1, col=c("violetred","black"))
out.hk <- scanone(fake.f2, chr=c(1,13), method="hk")
plot(out.hk, chr=c(1,13), add=TRUE, col="slateblue")
}

\seealso{ \code{\link[qtl]{scanone}},
  \code{\link[qtl]{summary.scanone}}, \code{\link[base]{par}},
  \code{\link[base]{colors}} }

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{hplot}
