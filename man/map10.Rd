\name{map10}
\alias{map10}
\non_function{}

\title{An example genetic map}

\description{ 
  A genetic map corresponding approximately to the mouse genome with
  a 10 cM marker spacing.
} 

\usage{data(map10)}

\format{
  An object of class \code{map}: a list whose components are vectors of
  marker locations.  This map approximates the mouse genome (see Rowe et
  al 1992), with 20 chromosomes (including the X 
  chromosome) and 148 markers at a 10 cM spacing.
}

\examples{
  data(map10)
  plot(map10)

  mycross <- sim.cross(map10, type="f2", n.ind=100)
}    

\seealso{ \code{\link[qtl]{sim.map}}, \code{\link[qtl]{plot.map}},
  \code{\link[qtl]{pull.map}} }

\keyword{datasets}
