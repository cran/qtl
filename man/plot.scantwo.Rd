\name{plot.scantwo}
\alias{plot.scantwo}

\title{Plot LOD scores for a two-dimensional genome scan}

\description{
  Plot the results of a two-dimensional, two-QTL genome scan.
}

\usage{
plot.scantwo(x, chr, incl.markers=FALSE, zlim,
             lower=c("cond-int","cond-add","joint"), nodiag=TRUE, 
             contours=FALSE, main, zscale=TRUE, \dots)
}			    

\arguments{
  \item{x}{An object of class \code{"scantwo"}, as output by
    \code{link[qtl]{scantwo}}.} 
  \item{chr}{Vector specifying which chromosomes to plot.}
  \item{incl.markers}{If FALSE, plot LOD scores on an evenly
    spaced grid (not including the results at the markers).}
  \item{zlim}{A vector of length 2 (optional), indicating the z limits
    for the joint LOD and the epistasis LOD scores, respectively.  If
    missing, the maximum limits are used for each.}
  \item{lower}{Indicates which LOD scores should be plotted in the lower
    triangle.  See the details below.}
  \item{nodiag}{If TRUE, suppress the plot of the scanone output
    (which is normally along the diagonal.)}
  \item{contours}{If TRUE, add contours to the plot, using a call to
    \code{\link[base]{contour}}.}
  \item{main}{An optional title for the plot.}
  \item{zscale}{If TRUE, a color scale is plotted at the right.}
  \item{\dots}{Ignored at this point.}
}

\value{None.}

\details{
  Uses \code{\link[base]{image}} to plot a grid of LOD scores.  The
  upper triangle contains the epistasis LOD scores; the content in the
  lower triangle depends on the input argument \code{lower}.
  The diagonal contains either all zeros or the main effects LOD scores
  (from \code{\link[qtl]{scanone}}).

  The \code{\link[qtl]{scantwo}} function calculates, for each pair of
  putative QTLs, \eqn{(q[1],q[2])}{(q1,q2)}, the likelihood undering the
  null model \eqn{L[0]}{L0}, the likelihood under each of the single-QTL
  models, \eqn{L(q[1])}{L(q1)} and \eqn{L(q[2])}{L(q2)}, the likelihood
  under an additive QTL model, \eqn{L[a](q[1],q[2])}{La(q1,q2)}, and the
  likelihood under a full QTL model (including QTL-QTL interaction),
  \eqn{L[f](q[1],q[2])}{Lf(q1,q2)}.  
  
  The upper triangle in plot from this function contains the epistasis
  LOD scores, \eqn{LOD[i] = log[10] L[f](q[1],q[2]) - log[10] L[a](q[1],q[2])}{%
    LODi = log10 Lf(q1,q2) - log10 Lf(q1,q2)}.

  If \code{lower="joint"}, the lower triangle will contain the joint LOD
  scores, \eqn{LOD[j] = log[10] L[f](q[1],q[2]) - log[10] L0}{%
    LODi = log10 Lf(q1,q2) - log10 L0}.

  If \code{lower="cond-int"}, the lower triangle will contain LOD scores
  comparing the two-QTL joint likelihood to the best of the single-QTL models,
  \eqn{LOD[ci] = log[10] L[f](q[1],q[2]) - log[10] max(L(q[1]), L(q[2]))}{%
    LODci = log10 Lf(q1,q2) - log10 max(L(q1), L(q2))}.

  If \code{lower="cond-add"}, the lower triangle will contain LOD scores
  comparing the two-QTL additive likelihood to the best of the single-QTL models,
  \eqn{LOD[ca] = log[10] L[a](q[1],q[2]) - log[10] max(L(q[1]), L(q[2]))}{%
    LODca = log10 La(q1,q2) - log10 max(L(q1), L(q2))}.  
    
  If a color scale is plotted (\code{zscale=TRUE}), the axis on the
  left indicates the scale for the epistasis LOD scores (the upper
  triangle), while the axis on the right indicates the scale for the
  lower triangle.  Note that the axis labels get screwed up
  if you change the size of the figure window; you'll need to redo the
  plot. 

  Note that the numbers plotted with the contours (if
  \code{contours=TRUE}) correspond to the scale for the lower triangle,
  and are likely incorrect for the upper triangle.  It appears to be
  difficult to fix this.   
}


\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2, step=5)
out.2dim <- scantwo(fake.f2, method="hk", chr=c(1,8,13,14))
plot(out.2dim)
plot(out.2dim, nodiag=FALSE)
plot(out.2dim, lower="joint", nodiag=FALSE)
}

\seealso{ \code{\link[qtl]{scantwo}},
  \code{\link[qtl]{summary.scantwo}}, \code{\link[qtl]{plot.scanone}} }

\author{Hao Wu, \email{hao@jax.org}; Karl W Broman,
  \email{kbroman@jhsph.edu} } 

\keyword{hplot}
