\name{lodint}
\alias{lodint}

\title{LOD support interval}

\description{
  Calculate a LOD support interval for a particular chromosomes.
}

\usage{lodint(results, chr, drop=1.5)}

\arguments{
 \item{results}{Output from \code{\link[qtl]{scanone}}.}
 \item{chr}{A chromosome ID.}
 \item{drop}{LOD units to drop to form the interval.}
}

\value{
  An object of class \code{scanone}, like the input, indicating the
  position with the maximum LOD, and indicating approximate endpoints
  for the LOD support interval.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)
listeria <- subset(listeria,chr=c(5,13))
listeria <- calc.genoprob(listeria, step=0.5)
out <- scanone(listeria, model="np")
lodint(out, 5, 1.5)
lodint(out, 13, 1.5)
}

\seealso{ \code{\link[qtl]{scanone}} }

\keyword{manip}
