\name{est.map}
\alias{est.map}

\title{Estimate genetic maps}

\description{
  Uses the Lander-Green algorithm (i.e., the hidden Markov model
  technology) to re-estimate the genetic map for an experimental cross.
}

\usage{
est.map(cross, error.prob=0, map.function=c("haldane","kosambi","c-f"),
        maxit=1000, tol=1e-5, sex.sp=TRUE, print.rf=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
 map function when converting genetic distances into recombination
 fractions.}
 \item{maxit}{Maximum number of EM iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps; this is 
 used only for the 4-way cross.}
 \item{print.rf}{Indicates whether to print initial and final estimates
   of the recombination fractions for each chromosome.  [This is mostly
   for debugging purposes.]}
}

\details{
  Calculations are done within the C function \code{est_map}.  The
  maximized log likelihood for each chromosome is saved as an
  attribute named \code{loglik}.  
}

\value{
  A \code{map} object; a list whose components (corresponding to
  chromosomes) are either vectors of marker positions (in cM) or
  matrices with two rows of sex-specific marker positions.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
newmap <- est.map(fake.f2)
logliks <- sapply(newmap,attr,"loglik")
plot.map(fake.f2,newmap)
fake.f2 <- replace.map(fake.f2,newmap)
}

\references{
  K Lange (1999) Numerical analysis for statisticians
  Springer-Verlag, New York. Sec 23.3.

  LR Rabiner (1989) A tutorial on hidden Markov models and selected
  applications in speech recognition.  Proceedings of the IEEE
  77:257-286.

  ES Lander, P Green (1987) Construction of multilocus genetic linkage
  maps in humans.  Proc Natl Acad Sci USA 84:2363-2367.
}  

\seealso{\code{\link[qtl]{plot.map}}, \code{\link[qtl]{est.rf}} }

\keyword{utilities}
