\name{pull.geno}
\alias{pull.geno}

\title{Pull out the genotype data from a cross}

\description{
  Pull out the genotype data from a cross object, as a single big matrix.
}

\usage{pull.geno(cross)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  A matrix of size n.ind x tot.mar.  The raw genotype data in the input
  cross object, with the chromosomes pasted together.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)
dat <- pull.geno(listeria)

# image of the genotype data
image(1:ncol(dat),1:nrow(dat),t(dat),ylab="Individuals",xlab="Markers",
      col=c("red","yellow","blue","green","violet"))
abline(v=cumsum(c(0,nmar(listeria)))+0.5)
abline(h=nrow(dat)+0.5)
}

\seealso{ \code{\link[qtl]{pull.map}} }

\keyword{manip}
