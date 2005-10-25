\name{est.map}
\alias{est.map}

\title{Estimate genetic maps}

\description{
  Uses the Lander-Green algorithm (i.e., the hidden Markov model
  technology) to re-estimate the genetic map for an experimental cross.
}

\usage{
est.map(cross, error.prob=0.0001,
        map.function=c("haldane","kosambi","c-f","morgan"), maxit=4000,
        tol=1e-4, sex.sp=TRUE, verbose=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
 \item{maxit}{Maximum number of EM iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps; this is 
 used only for the 4-way cross.}
 \item{verbose}{Logical; indicates whether to print initial and final
   estimates of the recombination fractions for each chromosome.}
}

\value{
  A \code{map} object; a list whose components (corresponding to
  chromosomes) are either vectors of marker positions (in cM) or
  matrices with two rows of sex-specific marker positions.
  The maximized log likelihood for each chromosome is saved as an
  attribute named \code{loglik}.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2,chr=18:19)}
newmap <- est.map(fake.f2)
logliks <- sapply(newmap, attr, "loglik")
plot.map(fake.f2, newmap)
fake.f2 <- replace.map(fake.f2, newmap)
}

\references{
  Lange, K. (1999) \emph{Numerical analysis for statisticians}.
  Springer-Verlag. Sec 23.3.

  Rabiner, L. R. (1989) A tutorial on hidden Markov models and selected
  applications in speech recognition.  \emph{Proceedings of the IEEE}
  \bold{77}, 257--286.

  Lander, E. S. and Green, P. (1987) Construction of multilocus genetic linkage
  maps in humans.  \emph{Proc. Natl. Acad. Sci. USA} \bold{84}, 2363--2367.
}  

\seealso{ \code{\link[qtl]{plot.map}}, \code{\link[qtl]{replace.map}},
  \code{\link[qtl]{est.rf}} }

\keyword{utilities}
