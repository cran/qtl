\name{summary.scantwo}
\alias{summary.scantwo}

\title{Summarize the results of a two-dimensional genome scan}

\description{
  Summarize the interesting aspects of the results of \code{scantwo}.
}

\usage{
summary.scantwo(object, thresholds, \dots)
}

\arguments{
 \item{object}{An object of class \code{scantwo}, the output of the
   function \code{\link[qtl]{scantwo}}.}
 \item{thresholds}{A vector of length three, giving LOD thresholds
   for the joint LOD, interaction LOD and single-QTL conditional LOD.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.scantwo}, to be printed by
  \code{print.summary.scantwo}.  Pairs of loci meeting the specified
  criteria are printed, with their joint LOD, interaction LOD, and the
  conditional LOD for each locus, along with single-point P-values
  calculated by the \eqn{\chi^2}{chi-square} approximation.  P-values
  are printed as \eqn{-\log_{10}(P)}{-log10(P)}. \cr

  If the input \code{scantwo} object does not include the results of
  \code{scanone}, the interaction and conditional LOD thresholds are
  ignored, and all pairs of loci for which the joint LOD exceeds its
  threshold are printed, though without their conditional LOD scores.
}

\details{
  For each pair of chromosomes, the pair of loci for which the joint LOD
  score is a maximum is considered.  The pair is printed only if its
  joint LOD score exceeds the joint threshold and either (a) the
  interaction LOD score exceeds its threshold or (b) both of the loci have
  conditional LOD scores that are above the conditional LOD threshold,
  where the conditional LOD score for locus \eqn{q_1}{q1}, LOD(\eqn{q_1
    | q_2}{q1 | q2}), is the \eqn{\log_{10}}{log10} likelihood ratio
  comparing the model with \eqn{q_1}{q1} and \eqn{q_2}{q2} acting
  additively to the model with \eqn{q_2}{q2} alone. 

  In the case the results of \code{\link[qtl]{scanone}} are not
  available, the maximum locus pair for each chromosome is printed
  whenever its joint LOD exceeds the joint LOD threshold.

  The criterion used in this summary is due to Gary Churchill and Saunak
  Sen, and deserves careful consideration and possible revision.
}

\author{Hao Wu, \email{hao@jax.org}; Karl W Broman,
  \email{kbroman@jhsph.edu} } 

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2, step=5)
out.2dim <- scantwo(fake.f2, method="hk", chr=c(1,6,13))
# All pairs of loci
summary(out.2dim)
# Pairs meeting specified criteria
summary(out.2dim, c(7, 3, 3))
# Pairs with both conditional LODs > 2
summary(out.2dim,c(0,1000,2))
# Pairs with interaction LOD is above 3
summary(out.2dim,c(0,3,1000))
}

\seealso{ \code{\link[qtl]{scantwo}}, \code{\link[qtl]{plot.scantwo}} }

\keyword{print}
