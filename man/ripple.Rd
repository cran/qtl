\name{ripple}
\alias{ripple}

\title{Compare marker orders}

\description{
  Investigate different marker orders for a given chromosome, comparing
  all possible permutations of a sliding window of markers.
}

\usage{
ripple(cross, chr, window=4, error.prob=0,
       map.function=c("haldane","kosambi","c-f"),
       maxit=4000, tol=1e-4, sex.sp=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{The chromosome to investigate.  Only one chromosome is allowed.}
 \item{window}{Number of markers to include in the sliding window of
   permuted markers.  Larger numbers result in the comparison of a
   greater number of marker orders, but will require a considerable
   increase in computation time.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
   map function when converting genetic distances into recombination
   fractions.}
 \item{maxit}{Maximum number of EM iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps; this is 
   used only for the 4-way cross.}
}

\details{
  Calculations are done by first constructing a matrix of marker orders
  and then making repeated calls to the R function
  \code{\link[qtl]{est.map}}.  Of course, it would be faster to do
  everything within C, but this was a lot easier to code.
}

\value{
  A matrix, given class \code{"ripple"}; the first set of columns are
  marker indices describing the order.  The last two columns are LOD
  scores (log base 10 likelihood ratios) comparing each order to the
  initial order and the estimated chromosome length for the given
  order.  Positive LOD scores indicate that the alternate order has more
  support than the original.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(badorder)
rip2 <- ripple(badorder, 1, 2)
summary(rip2)

badorder <- switch.order(badorder, 1, rip2[2,])
\dontrun{rip2p <- ripple(badorder, 1, 2)
summary(rip2p)}
}

\seealso{ \code{\link[qtl]{summary.ripple}}, \code{\link[qtl]{switch.order}},
  \code{\link[qtl]{est.map}}, \code{\link[qtl]{est.rf}} }

\keyword{utilities}
