\name{bayesint}
\alias{bayesint}

\title{Bayesian credible interval}

\description{
  Calculate an approximate Bayesian credible interval for a particular
  chromosome, using output from \code{\link[qtl]{scanone}}.
}

\usage{bayesint(results, chr, prob=0.95)}

\arguments{
 \item{results}{Output from \code{\link[qtl]{scanone}}.}
 \item{chr}{A chromosome ID.}
 \item{prob}{Probability coverage of the interval.}
}

\value{
  An object of class \code{scanone}, like the input, indicating the
  position with the maximum LOD, and indicating approximate endpoints
  for the Bayesian credible interval.
}

\details{
  We take \eqn{10^{LOD}}{10^LOD}, rescale it to have area 1, and then
  calculate the connected interval with density above some threshold
  and having coverage matching the target probability.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,4))}
hyper <- calc.genoprob(hyper, step=0.5)
out <- scanone(hyper, method="hk")
bayesint(out, chr=1)
bayesint(out, chr=4)
bayesint(out, chr=4, prob=0.99)
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{lodint}} }

\keyword{utilities}
