\name{plot.rf}
\alias{plot.rf}

\title{Plot recombination fractions}

\description{
  Plot a grid showing the recombination fractions for all pairs of
  markers, or of the LOD scores between pairs of markers.
}

\usage{
plot.rf(x, chr, which=c("both","lod","rf"), \dots)
}

\arguments{
  \item{x}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
  \item{chr}{Vector specifying which chromosomes to plot
  (optional)}
  \item{which}{Indicate whether to plot LOD scores, recombination
               fractions or both.}
  \item{\dots}{Ignored at this point.}
}

\value{None.}

\details{
  Uses \code{\link[base]{image}} to plot a grid showing the
  recombination fractions and/or LOD scores for all pairs of markers.
  (The LOD scores are for a test of \eqn{r = 1/2}.)
  If both are plotted, the recombination fractions are in the upper left
  triangle while the LOD scores are in the lower right triangle.
  Red corresponds to a large LOD or a small recombination fraction,
  while blue is the reverse.  Note that missing values appear in light
  gray. 

  Recombination fractions are transformed by \eqn{-4(\log_2{r}+1)}{%
    -4(log2(r)+1)} to make them on the same sort of scale as LOD
  scores.  Values of LOD or the transformed recombination fraction that
  are above 12 are set to 12.
}


\examples{
data(badorder)
badorder <- est.rf(badorder)
plot.rf(badorder)
}

\seealso{ \code{\link[qtl]{est.rf}}, \code{\link[base]{image}},
  \code{\link[qtl]{badorder}}, \code{\link[qtl]{ripple}} }

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{hplot}
