\name{top.errorlod}
\alias{top.errorlod}

\title{List genotypes with large error LOD scores}

\description{
  Prints those genotypes with error LOD scores above a specified
  cutoff.
}

\usage{
top.errorlod(cross, chr, cutoff=3, msg=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector specifying for which chromosomes the error LOD
   scores should be inspected.}
 \item{cutoff}{Only those genotypes with error LOD scores above this
   cutoff will be listed.}
 \item{msg}{If TRUE, print a message if there are no apparent errors.}
}

\value{
  A data.frame with 4 columns, whose rows correspond to the genotypes
  that are possibly in error.  The four columns give the chromosome
  number, individual number, marker name, and error LOD score.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,19))}
# run calc.genoprob, then calc.errorlod; use same error.prob
hyper <- calc.genoprob(hyper,error.prob=0.01)
hyper <- calc.errorlod(hyper,error.prob=0.01)

# plot the error LOD scores; print those above a specified cutoff 
plot.errorlod(hyper)
plot.errorlod(hyper,chr=1)
top.errorlod(hyper,cutoff=4)

# plot genotype data, flagging genotypes with error LOD > cutoff
plot.geno(hyper, chr=1, ind=160:200, cutoff=7, min.sep=2)

# Note that if markers are extremely tightly linked,
#     any recombination event will have a high error LOD
plot.geno(hyper, chr=1, ind=101:118, min.sep=2)
hyper$geno[[1]]$map
}

\seealso{ \code{\link[qtl]{calc.errorlod}},
  \code{\link[qtl]{plot.errorlod}} }

\keyword{print}
