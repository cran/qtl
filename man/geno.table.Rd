\name{geno.table}
\alias{geno.table}
\title{Create table of genotype distributions}
\description{
  Create table showing the observed numbers of individuals with each
  genotype at each marker.
}

\usage{
geno.table(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  A matrix of dimension [n.mar x n.gen] showing, for each marker, the number of
  individuals with each possible genotype, as well as the number that
  were not typed.
}

\examples{
data(listeria)
geno.table(listeria)
}

\seealso{\code{\link[qtl]{drop.markers}},
  \code{\link[qtl]{drop.nullmarkers}} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{manip}
