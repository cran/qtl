\name{plot.info}
\alias{plot.info}
\title{Plot the proportion of missing information in the genotype data}
\description{
  Plot a measure of the proportion of missing information in the
  genotype data.  
}

\usage{
plot.info(x,chr,which=c("both","entropy","variance"),return.result=FALSE,\dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Vector specifying which chromosomes to plot.}
 \item{which}{Indicates whether to plot the entropy version of the
   information, the variance version, or both.}
 \item{return.result}{If TRUE, returns the results (as a
   \code{data.frame}).}  
 \item{\dots}{Ignored at this point.}
}

\details{
  The entropy version of the missing information: for a single
  individual at a single genomic position, we measure the information as
  \eqn{H = \sum_g p_g \log p_g / \log n}, where \eqn{p_g} is the
  probability of the genotype \eqn{g}, and \eqn{n} is the number of
  possible genotypes.  This takes values between 0 and 1, assuming the
  value 0 when the genotypes are equally likely and 1 when one genotype
  has probability 0.  We calculate the missing information at a
  particular position as the average of \eqn{H} across individuals.  For
  an intercross, we don't scale by \eqn{\log n} but by the entropy in
  the case of genotype probabilities (1/4, 1/2, 1/4).
  
  The variance version of the missing information: we calculate the
  average, across individuals, of the variance of the genotype
  distribution at a particular locus, and scale by the maximum such
  variance. 

  Calculations are done in C, for the sake of speed in the presence of
  little thought about programming efficiency.
}  

\value{
  Calls \code{\link[qtl]{plot.scanone}} to plot the proportion of
  missing information in the genotype data.  If \code{return.result} is
  TRUE, the results are returned as a \code{data.frame} with class 
  \code{scanone}. 
}

\examples{
data(hyper)
hyper <- calc.genoprob(hyper,step=2,off.end=10)
plot.info(hyper)
}

\seealso{\code{\link[qtl]{plot.scanone}}}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
