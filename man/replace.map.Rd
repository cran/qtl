\name{replace.map}
\alias{replace.map}
\title{Replace the genetic map of a cross}

\description{
  Replace the map portion of a cross object.
}

\usage{replace.map(cross, map)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{map}{A list containing the new genetic map.  This must be the
   same length and with the same marker names as that contained in
   \code{cross}. }
}

\value{The input \code{cross} object with the genetic map replaced by
  the input \code{map}.}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2,fake.bc)
map <- pull.map(fake.f2)
fake.bc <- replace.map(fake.bc,map)
}

\seealso{\code{\link[qtl]{pull.map}} }

\keyword{manip}
