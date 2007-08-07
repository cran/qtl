\name{summary.scanoneboot}
\alias{summary.scanoneboot}

\title{Bootstrap confidence interval for QTL location}

\description{
  Calculates a bootstrap confidence interval for QTL location, using the
  bootstrap results from \code{\link[qtl]{scanoneboot}}.
}

\usage{
\method{summary}{scanoneboot}(object, prob=0.95, \dots)}
}

\arguments{
 \item{object}{Output from \code{\link[qtl]{scanoneboot}}.}
 \item{prob}{Desired coverage.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{scanone}, indicating the
  position with the maximum LOD, and indicating endpoints
  for the estimated bootstrap confidence interval.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
\dontrun{data(hyper)
hyper <- calc.genoprob(hyper, step=1, err=0.001)
bootoutput <- scanoneboot(hyper, chr=4)

hist(bootoutput, breaks=100)
summary(bootoutput)}
}

\seealso{ \code{\link[qtl]{scanoneboot}},
  \code{\link[qtl]{lodint}}, \code{\link[qtl]{bayesint}} }

\keyword{utilities}
