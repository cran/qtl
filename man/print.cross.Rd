\name{print.cross}
\alias{print.cross}


\title{Print a QTL experiment}

\description{
  Print method for a \code{cross} object.
}

\usage{
print.cross(x, \dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{\dots}{Ignored at this point.}
}

\value{
  A short warning message is printed, and then the function
  \code{\link[qtl]{summary.cross}} is called, and so a summary of the QTL
  object is printed, and the summary is returned.

  The main reason for including this function is that these
  \code{"cross"} objects are generally quite large and complex, and one
  often unintentionally types their name and gets the entire thing
  printed out.  The inclusion of this function makes typing the name of
  a \code{"cross"} object result in just the summary being printed.
}


\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(hyper)
hyper

print(hyper)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}}
          \code{\link[qtl]{summary.cross}} }

\keyword{print}
