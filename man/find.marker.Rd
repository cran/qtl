\name{find.marker}
\alias{find.marker}

\title{Find marker closest to a specified position}

\description{
  Find the genetic marker closest to a specified position on a
  chromosome. 
}

\usage{find.marker(cross, chr, pos)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector of chromosome identifiers, or a single such.}
 \item{pos}{A vector of cM positions.}
}

\value{
  A vector of marker names (of the same length as the input \code{pos}),
  corresponding to the markers nearest to the specified
  chromosomes/positions. 
}

\details{
  If the input \code{chr} has length one, it is expanded to the same
  length as the input \code{pos}.

  If multiple markers are exactly the same distance from the specified
  position, one is chosen at random from among those with the most
  genotype data.

  For a cross with sex-specific maps, the input positions are assumed to
  correspond to the female genetic map.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(listeria)
find.marker(listeria, 5, 28)
find.marker(listeria, c(1, 5, 13), c(81, 28, 26))
}

\seealso{ \code{\link[qtl]{find.flanking}}, \code{\link[qtl]{plot.pxg}},
  \code{\link[qtl]{find.pseudomarker}}, \code{\link[qtl]{effectplot}},
\code{\link[qtl]{find.markerpos}} }

\keyword{utilities}
