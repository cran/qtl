\name{fake.f2}
\alias{fake.f2}
\non_function{}
\title{Simulated data for an F2 intercross experiment}
\description{ 
  This is simulated data for an F2 intercross, obtained using
  \code{\link[qtl]{sim.cross}}.  The data were simulated using an additive
  model with three QTLs on chromosome 1 (at 30, 50 and 70 cM), and one QTL
  on chromosome 13 (at 30 cM).
} 
\usage{data(fake.f2)}
\format{
  The data is an object of class \code{cross}.  See
  \code{\link[qtl]{read.cross}} for details.
  There are 100 F2 individuals typed at 91 markers and with one phenotype. 
}
\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
}

\seealso{\code{\link[qtl]{sim.cross}}, \code{\link[qtl]{fake.bc}},
  \code{\link[qtl]{fake.4way}}, \code{\link[qtl]{listeria}},
  \code{\link[qtl]{hyper}} } 

\keyword{datasets}
