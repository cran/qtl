\name{clean}
\alias{clean}
\title{Remove intermediate calculations}
\description{
  Remove information on intermediate calculations from a cross object.
}

\usage{
clean(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  An object of class \code{cross}, with any intermediate calculations
  (such as is produced by \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{argmax.geno}} and \code{\link[qtl]{sim.geno}}) removed.
}

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2)
fake.f2 <- clean(fake.f2)
}

\seealso{\code{\link[qtl]{drop.qtlgeno}},
  \code{\link[qtl]{drop.nullmarkers}},
  \code{\link[qtl]{drop.markers}} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{manip}
