\name{pull.chr}
\alias{pull.chr}
\title{Pull out a subset of chromosomes from a cross}

\description{
  Pull out specified chromosomes from a cross object.
}

\usage{pull.chr(cross, chr)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector defining which chromosomes to pull.}
}

\value{The input \code{cross} object, but with the \code{geno} component
  containing only those chromosomes specified by the argument
  \code{chr}. }

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2.5n13 <- pull.chr(fake.f2,c("5","13"))
}

\seealso{\code{\link[qtl]{pull.map}}}

\keyword{manip}
