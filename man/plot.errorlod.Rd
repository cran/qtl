\name{plot.errorlod}
\alias{plot.errorlod}
\title{Plot grid of error LOD values}
\description{
  Plot a grid showing the LOD scores indicating which genotypes are
  likely to be in error, including lines at the chromosome boundaries.
}
\usage{
plot.errorlod(x,chr,\dots)
}
\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link[qtl]{read.cross}} for details.}
  \item{chr}{Indicates for which chromosomes the error LOD scores should
    be plotted.}
  \item{\dots}{Ignored at this point.}
}
\value{
  Uses \code{\link[base]{image}} to plot a grid with different shades
  of pixels to indicate which genotypes are most likely to be in
  error.  (Darker pixels have higher error LOD scores.) 
}

\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc,error.prob=0.01)
fake.bc <- calc.errorlod(fake.bc,error.prob=0.01)
plot.errorlod(fake.bc)
\testonly{data(fake.f2)}
\testonly{fake.f2 <- calc.errorlod(fake.f2)}
\testonly{plot.errorlod(fake.f2)}
}

\seealso{\code{\link[base]{image}}, \code{\link[qtl]{calc.errorlod}},
  \code{\link[qtl]{top.errorlod}}, \code{\link[qtl]{find.errors}},
  \code{\link[qtl]{plot.errors}} }

\references{SE Lincoln and ES Lander (1992) Systematic detection of
  errors in genetic linkage data.  Genomics 14:604-610}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
