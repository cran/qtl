\name{convert2ss}
\alias{convert2ss}
\title{Convert an intercross into a form to allow sex-specific maps}

\description{
  Convert a standard F2 intercross into a form which allows sex-specific
  maps (class \code{f2ss}).
}

\usage{convert2ss(cross)}

\arguments{
 \item{cross}{An object of class \code{cross} (must be an
   intercross). See \code{\link[qtl]{read.cross}} for details.} 
}

\value{The input \code{cross} object, but with the genetic map replaced
  by a pair (i.e., male and female) maps, and with the class
  changed to \code{f2ss}.}

\section{Warning}{
  The use of sex-specific maps for intercrosses is experimental.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)
listeria.ss <- convert2ss(listeria)
}

\seealso{ \code{\link[qtl]{est.map}}, \code{\link[qtl]{listeria}} }

\keyword{manip}
