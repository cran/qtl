\name{summary.ripple}
\alias{summary.ripple}

\title{Print summary of ripple results}

\description{
  Print marker orders, from the output of the function \code{ripple},
  for which the log10 likelihood relative to the initial order is above
  a specified cutoff.
}

\usage{
summary.ripple(object, lod.cutoff = -1, \dots)
}

\arguments{
 \item{object}{An object of class \code{ripple}, the output of the
   function \code{\link[qtl]{ripple}}.}
 \item{lod.cutoff}{Only marker orders with LOD score (relative to the
   initial order) above this cutoff will be displayed.} 
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.ripple}, whose rows correspond to
  marker orders with likelihood within some cutoff of the initial order.
  If no marker order, other than the initial one, has likelihood within
  the specified range, the initial and next-best orders are returned.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{ % I don't run the following, since they are run in ripple.Rd
\dontrun{data(badorder)
rip2 <- ripple(badorder, 1, 2)
summary(rip2)

badorder <- switch.order(badorder, 1, rip2[2,])
rip2p <- ripple(badorder, 1, 2)
summary(rip2p)}
}

\seealso{ \code{\link[qtl]{ripple}}, \code{\link[qtl]{est.map}},
  \code{\link[qtl]{est.rf}} }

\keyword{print}
