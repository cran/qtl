\name{summary.ripple}
\alias{summary.ripple}

\title{Print summary of ripple results}

\description{
  Print marker orders, from the output of the function \code{ripple},
  for which the likelihood is within some cutoff of the initial order.
}

\usage{
summary.ripple(object, lod.cutoff=1, \dots)
}
\arguments{
 \item{object}{An object of class \code{ripple}, the output of the
   function \code{\link[qtl]{ripple}}.}
 \item{lod.cutoff}{Only marker orders with LOD score (relative to the
   initial order) above the negative of this cutoff will be displayed.} 
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.ripple}, whose rows correspond to
  marker orders with likelihood within some cutoff of the initial order.
  If no marker order, other than the initial one, has likelihood within
  the specified range, the initial and next-best orders are returned.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
\dontrun{data(fake.f2)}
\dontrun{output <- ripple(fake.f2,7,2)}
\dontrun{summary(output,2)}
}

\seealso{\code{\link[qtl]{ripple}}, \code{\link[qtl]{est.map}},
  \code{\link[qtl]{est.rf}} }

\keyword{utilities}
