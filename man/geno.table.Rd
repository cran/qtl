\name{geno.table}
\alias{geno.table}

\title{Create table of genotype distributions}

\description{
  Create table showing the observed numbers of individuals with each
  genotype at each marker, including P-values from chi-square tests
  for Mendelian segregation.
}

\usage{
geno.table(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  A matrix containing, for each marker, the number of individuals with
  each possible genotype, as well as the number that were not typed.
  The first column gives the chromosome ID, and the last column gives
  P-values from chi-square tests of Mendelian segregation.  In the case
  of partially informative genotypes, P-values are not provided.
}

\examples{
data(listeria)
geno.table(listeria)
}

\seealso{ \code{\link[qtl]{summary.cross}},
  \code{\link[qtl]{drop.markers}}, \code{\link[qtl]{drop.nullmarkers}} }

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{utilities}
