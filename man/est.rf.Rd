\name{est.rf}
\alias{est.rf}

\title{Estimate pairwise recombination fractions}

\description{
  Estimate the sex-averaged recombination fraction between all pairs of
  typed markers.
}

\usage{
est.rf(cross, maxit=1000, tol=1e-6)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{maxit}{Maximum number of iterations for the EM algorithm (not
   used with backcrosses).}
 \item{tol}{Tolerance for determining convergence (not used with
   backcrosses).} 
}

\value{
  The \code{cross} object in the input is returned with a component,
  \code{rf}, added.  This is a matrix of size (tot.mar x tot.mar).  The
  diagonal contains the number of typed meioses per marker, the lower
  triangle contains the estimated recombination fractions, and the upper
  triangle contains the LOD scores (testing rf = 0.5).
}

\details{
  For a backcross, one can simply count recombination events.  For an
  intercross or 4-way cross, a version of the EM algorithm must be used
  to estimate recombination fractions. (Since, for example, in an
  intercross individual that is heterozygous at two loci, it is not
  known whether it shows 0 or 2 recombinations.)  Note that, for the
  4-way cross, we estimate sex-averaged recombination fractions.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2 <- est.rf(fake.f2)
plot.rf(fake.f2)
}

\seealso{\code{\link[qtl]{plot.rf}}, \code{\link[qtl]{est.map}} }

\keyword{utilities}
