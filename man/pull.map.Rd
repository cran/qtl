\name{pull.map}
\alias{pull.map}
\title{Pull out the genetic map from a cross}

\description{
  Pull out the map portion of a cross object.
}

\usage{pull.map(cross)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{ A list containing the genetic map.}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
map <- pull.map(fake.f2)
}

\seealso{\code{\link[qtl]{replace.map}},
  \code{\link[qtl]{pull.chr}} }

\keyword{manip}
