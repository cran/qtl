\name{drop.markers}
\alias{drop.markers}
\title{Drop a set of markers}
\description{
  Drop a vector of markers the data matrices and genetic maps.
}

\usage{
drop.markers(cross,markers)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{markers}{A character vector of marker names.}
}

\value{
  The same object in the input, with any markers in the vector
  \code{markers} removed from the genotype data matrices, genetic maps,
  and, if applicable, the results of \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{argmax.geno}}, \code{\link[qtl]{sim.geno}},
  and \code{\link[qtl]{est.rf}}. 
}

\examples{
data(listeria)
listeria2 <- drop.markers(listeria,c("D10M44","D1M3","D1M75","D1M155"))
}

\seealso{\code{\link[qtl]{drop.nullmarkers}}, \code{\link[qtl]{geno.table}},
  \code{\link[qtl]{clean}}, \code{\link[qtl]{drop.qtlgeno}} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{manip}
