\name{fake.bc}
\alias{fake.bc}
\non_function{}
\title{Simulated data for a backcross experiment}
\description{ 
  This is simulated data for a backcross, obtained using
  \code{\link[qtl]{sim.cross}}.  There are two correlated phenotypes, each
  due to an additive model with three QTLs on chromosome 1 (at 30, 50 and
  70 cM), and one QTL on chromosome 9 (at 30 cM).  The second phenotype
  is also influenced by a QTL on chromosome 11 (at -5 cM) and two QTLs
  on chromosome 13 (at 10 and 30 cM).
} 
\usage{data(fake.bc)}
\format{
  The data is an object of class \code{cross}.  See
  \code{\link[qtl]{read.cross}} for details.
  There are 400 backcross individuals typed at 91 markers and with two
  phenotypes.
}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.bc)
}

\seealso{\code{\link[qtl]{sim.cross}}, \code{\link[qtl]{fake.4way}},
  \code{\link[qtl]{fake.f2}}, \code{\link[qtl]{listeria}},
  \code{\link[qtl]{hyper}} }

\keyword{datasets}
