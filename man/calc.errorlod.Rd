\name{calc.errorlod}
\alias{calc.errorlod}
\title{Identify likely genotyping errors}
\description{
  Calculates a LOD score for each genotype, indicating which ones are
  likely to be in error.  
}
\usage{
calc.errorlod(cross, error.prob=0.01, map.function=c("haldane","kosambi","c-f"))
}
\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype)}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
 map function when converting genetic distances into recombination
 fractions.}
}
\details{
  Uses the function \code{\link[qtl]{calc.genoprob}} to calculate conditional
  genotype probabilities given observed marker genotypes, and then
  calculates the LOD score described by Lincoln and Lander (1992) for
  each genotype, indicating likely errors.  The actual calculations are
  done in the C program \code{calc_errorlod}.

  An attribute "error.prob" is set to the value of the corresponding
  argument, for later reference. 
}
\value{
  The \code{cross} object in the input is returned with a component,
  \code{errorlod}, added to each component of \code{cross$geno}.  The
  \code{errorlod} component is a matrix of size (n.ind x n.mar).
}
\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\references{SE Lincoln and ES Lander (1992) Systematic detection of
  errors in genetic linkage data.  Genomics 14:604-610}

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2,error.prob=0.01)
fake.f2 <- calc.errorlod(fake.f2,error.prob=0.01)
plot.errorlod(fake.f2)
\testonly{data(fake.bc)}
\testonly{fake.bc <- calc.errorlod(fake.bc)}
\testonly{plot.errorlod(fake.bc)}
}

\seealso{\code{\link[qtl]{plot.errorlod}},
  \code{\link[qtl]{top.errorlod}},
  \code{\link[qtl]{find.errors}} }

\keyword{utilities}
