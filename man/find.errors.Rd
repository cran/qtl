\name{find.errors}
\alias{find.errors}

\title{Identify likely genotyping errors}

\description{
  Compares the observed genotypes to the most likely sequence of
  underlying genotypes given the observed marker data (calculated in
  \code{\link[qtl]{argmax.geno}}), to identify genotypes likely to be in
  error. 
}

\usage{
find.errors(cross, chr, error.prob=0.01,
            map.function=c("haldane","kosambi","c-f"),msg=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector specifying for which chromosomes the calculations
   should be performed.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).} 
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
 map function when converting genetic distances into recombination
 fractions.} 
 \item{msg}{If TRUE, print a message if there are no apparent errors.}
}

\details{
  The most likely sequence of underlying genotypes is calculated in
  \code{\link[qtl]{argmax.geno}}.  This function simply pulls out those
  genotypes where the observed genotype differs from that appearing in
  the \code{\link[qtl]{argmax.geno}} results.
}

\value{
  A matrix with 5 columns, whose rows correspond to the genotypes that
  are possibly in error.  The five columns give the chromosome number,
  individual number, marker name, observed genotype and most likely
  genotype.  
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
output <- find.errors(fake.f2)
}

\seealso{\code{\link[qtl]{argmax.geno}},
         \code{\link[qtl]{plot.errors}}}

\keyword{utilities}
