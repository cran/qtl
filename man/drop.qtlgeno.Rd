\name{drop.qtlgeno}
\alias{drop.qtlgeno}
\title{Remove QTL genotypes}
\description{
  Remove QTLs from the genotype data and genetic maps of a simulated QTL
  cross object.
}

\usage{
drop.qtlgeno(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  An object of class \code{cross}, with any QTLs (markers with names
  beginning "QTL") removed from the genotype data and genetic maps.
}

\examples{
data(fake.f2)
fake <- sim.cross(pull.map(fake.f2),keep.qtlgeno=TRUE,model=c(1,5,0.5,0.6))
fake2 <- drop.qtlgeno(fake)
}

\seealso{\code{\link[qtl]{sim.cross}}, \code{\link[qtl]{clean}},
  \code{\link[qtl]{drop.nullmarkers}},
  \code{\link[qtl]{drop.markers}} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{manip}
