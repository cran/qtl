\name{drop.nullmarkers}
\alias{drop.nullmarkers}
\title{Drop markers without any genotype data}
\description{
  Drop markers, from the data matrices and genetic maps, that have no
  genotype data.
} 

\usage{
drop.nullmarkers(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  The same object in the input, with any markers lacking genotype
  data removed from the genotype data matrices, genetic maps, and, if
  applicable, the results of \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{argmax.geno}}, \code{\link[qtl]{sim.geno}},
  and \code{\link[qtl]{est.rf}}.
}

\examples{
data(listeria)
listeria2 <- drop.nullmarkers(listeria)
}

\seealso{\code{\link[qtl]{drop.markers}}, \code{\link[qtl]{clean}},
  \code{\link[qtl]{drop.qtlgeno}}, \code{\link[qtl]{geno.table}} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{manip}
