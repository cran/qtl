\name{switch.order}
\alias{switch.order}
\title{Switch the order of markers on a chromosome}

\description{
  Switch the order of markers on a specified chromosome to a specified
  new order.
}

\usage{switch.order(cross, chr, order)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{The chromosome for which the marker order is to be
   switched.}
 \item{order}{A vector of indices defining the new marker order. The
   vector may have length two more than the number of markers, for easy
   in use with the output of the function \code{\link[qtl]{ripple}}.}
}

\value{
  The input \code{cross} object, but with the marker order on the
  specified chromosome updated.

  The output object is "cleaned" of any intermediate calculations (such
  as from \code{\link[qtl]{calc.genoprob}}, because these values have no
  meaning with a change in marker order.

  We re-estimate the genetic map for the relevant chromosome.
}


\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2 <- switch.order(fake.f2,1,c(1,3,2,4:7))
}

\seealso{\code{\link[qtl]{ripple}}, \code{\link[qtl]{clean}}}

\keyword{manip}
