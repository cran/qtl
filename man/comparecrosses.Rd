\name{comparecrosses}
\alias{comparecrosses}
\title{Compare two cross objects.}

\description{
  Verify that two objects of class \code{cross} have identical classes,
  chromosomes, markers, genotypes, genetic maps, and phenotypes.
}

\usage{comparecrosses(cross1, cross2)}

\arguments{
 \item{cross1}{An object of class \code{cross} (must be an
   intercross). See \code{\link[qtl]{read.cross}} for details.} 
 \item{cross2}{An object of class \code{cross} (must be an
   intercross). See \code{\link[qtl]{read.cross}} for details.} 
}

\value{None.}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)
comparecrosses(listeria, listeria)
}

\seealso{ \code{\link[qtl]{summary.cross}} }

\keyword{manip}
