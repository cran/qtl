\name{c.cross}
\alias{c.cross}

\title{Combine data for QTL experiments}

\description{
  Concatenate the data for multiple QTL experiments.
}

\usage{c.cross(\dots)}

\arguments{
 \item{\dots}{A set of objects of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.  These must all be of the
   same cross type, with the same phenotypes and the same chromosomes
   and markers.}
}

\value{The concatenated input, as a \code{cross} object.}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(fake.f2)
junk <- fake.f2
junk <- c(fake.f2,junk)
}

\seealso{ \code{\link[qtl]{subset.cross}} }

\keyword{manip}
