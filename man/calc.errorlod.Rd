\name{calc.errorlod}
\alias{calc.errorlod}

\title{Identify likely genotyping errors}

\description{
  Calculates a LOD score for each genotype, measuring the evidence for
  genotyping errors. 
}

\usage{
calc.errorlod(cross, error.prob=0.01,
              map.function=c("haldane","kosambi","c-f"))
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype)}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
 map function when converting genetic distances into recombination
 fractions.}
}

\details{
  Uses the results of the function \code{\link[qtl]{calc.genoprob}}
  (conditional genotype probabilities given multipoint marker data), and
  then calculates, for each individual at each marker, a LOD score
  measuring the strength of evidence for a genotyping error, as
  described by Lincoln and Lander (1992). 

  If results of \code{\link[qtl]{calc.genoprob}} are not contained in
  the input \code{cross} object, it is called with arguments \code{step
  = off.end = 0}.  If it had been run, but with a different value of
  \code{error.prob} then specified in the argument, it is called again
  with the specified value, and using the available \code{"step"} and
  \code{"off.end"} attributes, and the only genotype probabilities are
  replaced.

  Note that values below 3 are generally not interesting.  Also note
  that if markers are extremely tightly linked, \emph{recombination
  events} can give large error LOD scores.  The error LOD scores should
  not be trusted blindly, but should be viewed as a tool for identifying
  genotypes to study further.  

  Use \code{\link[qtl]{top.errorlod}} to print all genotypes with error
  LOD scores above a specified threshold,
  \code{\link[qtl]{plot.errorlod}} to plot the error LOD scores for
  specified chromosomes, and \code{\link[qtl]{plot.geno}} to view the
  observed genotype data with likely errors flagged.
}

\value{
  The input \code{cross} object is returned with a component,
  \code{errorlod}, added to each component of \code{cross$geno}.  The
  \code{errorlod} component is a matrix of size (n.ind x n.mar).  An
  attribute \code{"error.prob"} is set to the value of the corresponding 
  argument, for later reference.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\references{
  Lincoln, S. E. and Lander, E. S. (1992) Systematic detection of
  errors in genetic linkage data.  \emph{Genomics} \bold{14}, 604--610.
}

\examples{
data(hyper)
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

\seealso{ \code{\link[qtl]{plot.errorlod}},
  \code{\link[qtl]{top.errorlod}} } 


\keyword{utilities}
