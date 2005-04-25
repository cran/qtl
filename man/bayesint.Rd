\name{bayesint}
\alias{bayesint}

\title{Bayesian probability interval}

\description{
  Calculate a Bayesian probability interval for a particular chromosome,
  using output from scanone.
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
  for the Bayesian probability interval.
}

\details{
  We take \eqn{10^{LOD}}{10^LOD}, rescale it to have area approximately
  1, and the look for an interval whose area is the target probability.
  
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)
listeria <- subset(listeria,chr=c(5,13))
listeria <- calc.genoprob(listeria, step=0.5)
out <- scanone(listeria, model="np")
bayesint(out, 5)
bayesint(out, 13)
bayesint(out, 13, 0.99)
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{lodint}} }

\keyword{manip}
