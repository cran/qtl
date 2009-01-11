\name{summary.cross}
\alias{summary.cross}

\title{Print summary of QTL experiment}

\description{
  Print summary information about a \code{cross} object.
}

\usage{
\method{summary}{cross}(object, \dots)
}

\arguments{
 \item{object}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.cross}
  containing a variety of summary information about the cross (this is
  generally printed automatically).
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
summary(fake.f2)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}},
  \code{\link[qtl]{nind}},
  \code{\link[qtl]{nmar}},
  \code{\link[qtl]{nchr}},
  \code{\link[qtl]{totmar}},
  \code{\link[qtl]{nphe}} }

\keyword{print}
