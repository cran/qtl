\name{convert2ss}
\alias{convert2ss}
\title{Convert an intercross into a form to allow sex-specific maps}

\description{
  Convert a standard F2 intercross (class "f2") into a form which allows
  sex-specific maps (class "f2ss").
}

\usage{convert2ss(cross)}

\arguments{
 \item{cross}{An object of class \code{cross}, must be an intercross.
   See \code{\link[qtl]{read.cross}} for details.}
}

\value{The input \code{cross} object, but with the genetic map replaced
  by a pair (i.e., male and female) maps, and with the "f2" class
  changed to "f2ss."}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2ss <- convert2ss(fake.f2)
}

\seealso{\code{\link[qtl]{est.map}},\code{\link[qtl]{fake.f2}}}

\keyword{manip}
