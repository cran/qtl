\name{plot.rf}
\alias{plot.rf}
\title{Plot recombination fractions}
\description{
  Plot a grid showing the recombination fractions for all pairs of
  markers, or of the LOD scores between pairs of markers.
}
\usage{
plot.rf(x, chr, which="both", \dots)
}
\arguments{
  \item{x}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
  \item{chr}{Vector specifying which chromosomes to plot
  (optional)}
  \item{which}{Indicate whether to plot LOD scores, recombination
               fractions or both (\code{lod}, \code{rf} or 
               \code{both}).}
  \item{\dots}{Ignored at this point.}
}
\value{
  Uses \code{\link[base]{image}} to plot a grid showing either the
  recombination fractions or LOD scores for all pairs of markers or
  both (in which case LOD scores are in the lower triangle and
  recombination fractions are in the upper triangle).
  Note that missing values appear in gray.
}

\examples{
data(fake.bc)
fake.bc <- est.rf(fake.bc)
plot.rf(fake.bc)
plot.rf(fake.bc,which="lod")
plot.rf(fake.bc,c(1,5),which="rf")
}

\seealso{\code{\link[base]{image}}, \code{\link[qtl]{est.rf}} }

	 

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
