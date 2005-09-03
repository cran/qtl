\name{arithscan}
\alias{+.scanone}
\alias{-.scanone}

\title{Arithmetic Operators for genome scan results}

\description{
  These binary operators perform arithmetic on genome scan result
  objects.
}

\usage{
\special{scan1+scan2}
\special{scan1-scan2}
}

\arguments{
  \item{scan1, scan2}{Genome scan results on the same set of chromosomes
    and markers.}
}

\details{
  This is used to calculate the sum or difference of LOD scores of two
  genome scan
  results. The basic arithmetic operators were overrided to work on
  objects of class \code{scanone}, \code{scantwo}, and \code{scanqtl}.

  It currently only works for \code{scanone} results.
}

\value{
  The same data structure as the input objects, with LOD scores
  recalcuated.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu}}

\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc)
# genome scans without covariates
out.nocovar <- scanone(fake.bc)
# genome scans with covariates
ac <- fake.bc$pheno[,c("sex","age")]
ic <- fake.bc$pheno[,"sex"]
out.covar <- scanone(fake.bc, pheno.col=1,
                     addcovar=ac, intcovar=ic)
# plot the difference of with and without additive covariate
plot(out.covar-out.nocovar)
}

\keyword{arith}
