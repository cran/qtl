\name{clean}
\alias{clean}

\title{Remove derived data}

\description{
  Remove any intermediate calculations from a cross object.
}

\usage{
clean(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  The input object, with any intermediate calculations
  (such as is produced by \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{argmax.geno}} and \code{\link[qtl]{sim.geno}})
  removed.
}

\examples{
data(fake.f2)
names(fake.f2$geno)
fake.f2 <- calc.genoprob(fake.f2)
names(fake.f2$geno)
fake.f2 <- clean(fake.f2)
names(fake.f2$geno)
}

\seealso{ \code{\link[qtl]{drop.nullmarkers}},
  \code{\link[qtl]{drop.markers}} }

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{manip}
