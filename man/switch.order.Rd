\name{switch.order}
\alias{switch.order}
\title{Switch the order of markers on a chromosome}

\description{
  Switch the order of markers on a specified chromosome to a specified
  new order.
}

\usage{switch.order(cross, chr, order, error.prob=0,
                    map.function=c("haldane","kosambi","c-f","morgan"))}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{The chromosome for which the marker order is to be
   switched.}
 \item{order}{A vector of numeric indices defining the new marker
   order. The vector may have length two more than the number of
   markers, for easy in use with the output of the function
   \code{\link[qtl]{ripple}}.} 
 \item{error.prob}{Assumed genotyping error rate (passed to
   \code{\link[qtl]{est.map}}).} 
 \item{map.function}{Map function to be used (passed to
   \code{\link[qtl]{est.map}}).} 
}

\value{
  The input \code{cross} object, but with the marker order on the
  specified chromosome updated, and with any derived data removed
  (except for recombination fractions, if present, which are not
  removed);  the genetic map for the relevant chromosome is
  re-estimated. 
}


\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(fake.f2)
fake.f2 <- switch.order(fake.f2, 1, c(1,3,2,4:7))
}

\seealso{ \code{\link[qtl]{ripple}}, \code{\link[qtl]{clean}} }

\keyword{manip}
