\name{plot.missing}
\alias{plot.missing}
\title{Plot grid of missing genotypes}

\description{
  Plot a grid showing which genotypes are missing, including lines at
  the chromosome boundaries.  
}

\usage{
plot.missing(x,chr,reorder=F,main="Missing genotypes",\dots)
}			    

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link[qtl]{read.cross}} for details. }
  \item{chr}{Vector of chromosomes to plot.  If missing, plot all
    chromosomes.}
  \item{reorder}{Specify whether to reorder individuals according to
  their phenotypes. 
  \tabular{ll}{
    FALSE \tab Don't reorder \cr
    TRUE  \tab Reorder according to the sum of the phenotypes \cr
    n \tab Reorder according to phenotype n}}
  \item{main}{Title to place on plot.}
  \item{\dots}{Ignored at this point.}
}

\value{
  Uses \code{\link[base]{image}} to plot a grid with black pixels where the
  genotypes are missing.  For intercross and 4way cross data, pink
  pixels are plotted for the partially missing genotypes.
}

\examples{
data(fake.f2)
plot.missing(fake.f2)
\testonly{data(fake.bc)}
\testonly{plot.missing(fake.bc)}
}

\seealso{\code{\link[base]{image}}, \code{\link[qtl]{plot.cross}} }


\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
