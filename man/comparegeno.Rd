\name{comparegeno}
\alias{comparegeno}

\title{Compare individuals' genotype data}

\description{
  Count proportion of matching genotypes between all pairs of
  individuals, to look for unusually closely related individuals.
}

\usage{comparegeno(cross, what=c("proportion","number"))}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{what}{Indicates whether to return the proportion or number of matching
genotypes.}
}

\value{
  A matrix whose (i,j)th element is the proportion or number of matching
genotypes for individuals i and j.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)

output <- comparegeno(listeria)

n.ind <- nind(listeria)
image(1:n.ind,1:n.ind,output,col=gray((0:99)/99),breaks=seq(0,1,len=101))
}

\seealso{ \code{\link[qtl]{nmissing}} }

\keyword{manip}
