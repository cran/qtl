\name{plot.cross}
\alias{plot.cross}
\title{Plot various features of a cross object}
\description{
  Plot grid of the missing genotypes, genetic map, and histograms of
  phenotypes for the data from an experimental cross.
}

\usage{
plot.cross(x,\dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{\dots}{Ignored at this point.}
}

\value{
  Calls \code{\link[qtl]{plot.missing}} and \code{\link[qtl]{plot.map}} and
  \code{\link[base]{hist}} to plot the missing genotypes, genetic map, and
  histograms of all phenotypes.
}

\examples{
data(fake.bc)
plot.cross(fake.bc)
plot(fake.bc)
}

\seealso{\code{\link[qtl]{plot.missing}}, \code{\link[qtl]{plot.map}}, 
         \code{\link[base]{hist}}}			   

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
