\name{fake.4way}
\alias{fake.4way}
\non_function{}
\title{Simulated data for a 4-way cross experiment}
\description{ 
  This is simulated data for a 4-way cross, obtained using
  \code{\link[qtl]{sim.cross}}.  There are two phenotypes (including sex, for
  which 0=female and 1=male).  The quantitative phenotype is affected by
  three QTLs: two on chromosome 2 at positions 10 and 25 cM on the female
  genetic map, and one on chromosome 7 at position 40 cM on the female map.
} 
\usage{data(fake.4way)}
\format{
  The data is an object of class \code{cross}.  See
  \code{\link[qtl]{read.cross}} for details.
  There are 250 individuals typed at 157 markers, including 8 on the X chromosome.
}
\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.4way)
}

\seealso{\code{\link[qtl]{sim.cross}}, \code{\link[qtl]{fake.bc}},
  \code{\link[qtl]{fake.f2}}, \code{\link[qtl]{listeria}},
  \code{\link[qtl]{hyper}} }

\keyword{datasets}
