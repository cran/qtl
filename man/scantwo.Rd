\name{scantwo}
\alias{scantwo}

\title{Two-dimensional genome scan with a two-QTL model}

\description{
  Perform a two-dimensional genome scan with a two-QTL model, with
  possible allowance for covariates. \bold{Warning}: two-dimensional
  genome scans on a dense grid can take a great deal of computer time
  and memory. 
}

\usage{
scantwo(cross, chr, pheno.col=1,
        method=c("em","imp","hk","mr","mr-imp","mr-argmax"), 
        addcovar=NULL, intcovar=NULL, weights=NULL,
        run.scanone=TRUE, incl.markers=FALSE, maxit=4000, tol=1e-4,
        trace=TRUE, n.perm)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Vector indicating the chromosomes for which LOD scores
   should be calculated.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.}
 \item{method}{Indicates whether to use standard interval mapping (ie
   the EM algorithm), imputation, Haley-Knott regression, or marker
   regression.  Marker regression is performed either by dropping
   individuals with missing genotypes (\code{"mr"}), or by first filling
   in missing data using a single imputation (\code{"mr-imp"}) or by the
   Viterbi algorithm (\code{"mr-argmax"}).}
 \item{addcovar}{Additive covariates.}
 \item{intcovar}{Interactive covariates (interact with QTL genotype).}
 \item{weights}{Optional weights of individuals.  Should be either NULL
   or a vector of length n.ind containing positive weights.}
 \item{run.scanone}{If TRUE, run the function \code{\link[qtl]{scanone}}
   and place the results on the diagonal.}
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.}
 \item{maxit}{Maximum number of iterations in the EM algorithm; used
   only with method \code{"em"}.}
 \item{tol}{Tolerance value for determining convergence in the EM
   algorithm; used only with method \code{"em"}.}
 \item{trace}{If TRUE, display information about the progress of
   calculations.  For method \code{"em"}, if \code{trace} is an integer
   above 1, further details on the progress of the algorithm will be
   displayed.}
 \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
}

\details{
  The interval mapping (i.e. EM algorithm) and Haley-Knott regression
  methods require that multipoint genotype probabilities are first
  calculated using \code{\link[qtl]{calc.genoprob}}.  The imputation
  method uses the results of \code{\link[qtl]{sim.geno}}. 

  The method \code{em} is standard interval mapping by the EM algorithm
  (Dempster et al. 1977; Lander and Botstein 1989).
  Marker regression is simply linear regression of phenotypes on marker
  genotypes  (individuals with missing genotypes are discarded).
  Haley-Knott regression uses the regression of phenotypes on multipoint
  genotype probabilities.  The imputation method uses the pseudomarker
  algorithm described by Sen and Churchill (2001).
  
  Individuals with missing phenotypes are dropped.

  In the presence of covariates, the full model is 
  \deqn{y = \mu + \beta_{q1} + \beta_{q2} + \beta_{q1 \times q2} +
    A \gamma + Z \delta_{q1} + Z \delta_{q2} + Z \delta_{q1 \times q2} +
    \epsilon}{y = m + b[q1] + b[q2] + b[q1 x q2] + A g + Z d[q1] + Z
    d[q2] + Z d[q1 x q2] + e}
  where \emph{q1} and \emph{q2} are the unknown QTL genotypes at two
  locations, \emph{A} is a matrix of covariates, and \emph{Z} is a
  matrix of covariates that interact with QTL genotypes.  The columns of
  \emph{Z} are forced to be contained in the matrix \emph{A}.

  We calculate LOD scores testing comparing the full model to each of
  two alternatives.  The joint LOD score compares the full model to the
  following null model:
  \deqn{y = \mu + A \gamma + \epsilon}{y = m + A g + e}
  The epistasis LOD score compares the full model to the following
  additive model: 
  \deqn{y = \mu + \beta_{q1} + \beta_{q2} + A \gamma + Z \delta_{q1} + Z
    \delta_{q2} + \epsilon}{y = m + b[q1] + b[q2] + A g + Z d[q1] + Z
    d[q2] + e}

  In the case that \code{n.perm} is specified, the R function
  \code{scantwo} is called repeatedly.
}

\section{X chromosome}{
The X chromosome must be treated specially in QTL mapping.

As in \code{\link[qtl]{scanone}}, if both males and females are
included, male hemizygotes are allowed to be different from female
homozygotes, and the null hypothesis must be changed in order to ensure
that sex- or pgm-differences in the phenotype do not results in spurious
linkage to the X chromosome.  
}

\value{
  If \code{n.perm} is missing, the function returns a list with 
  class \code{"scantwo"} and containing three components.  The first
  component is a matrix of dimension [tot.pos x tot.pos] whose upper
  triangle contains the epistasis LOD scores and whose lower triangle
  contains the joint LOD scores.  If \code{run.scanone=TRUE}, the
  diagonal contains the results of \code{\link[qtl]{scanone}}. The
  second component of the output is a data.frame indicating the
  locations at which the two-QTL LOD scores were calculated.  The first
  column is the chromosome identifier, the second column is the position
  in cM, the third column is a 1/0 indicator for ease in later
  pulling out only the equally spaced positions, and the fourth column
  indicates whether the position is on the X chromosome or not.  The
  final component is a version of the results of
  \code{\link[qtl]{scanone}} including sex and/or cross direction as
  additive covariates, which is needed for a proper calculation of
  conditional LOD scores.

  If \code{n.perm} is specified, the function returns a matrix with two
  columns, containing the maximum joint and epistasis LOD scores, across
  a two-dimensional scan, for each of the permutation replicates.
}

\references{
  Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for
  quantitative trait mapping.  \emph{Genetics} \bold{138}, 963--971.

  Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum
  likelihood from incomplete data via the EM algorithm.  \emph{J. Roy.
    Statist. Soc.} B, \bold{39}, 1--38.

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Lander, E. S. and Botstein, D. (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps.  \emph{Genetics}
  \bold{121}, 185--199.

  Sen, S. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.

  Soller, M., Brody, T. and Genizi, A. (1976) On the power of experimental
  designs for the detection of linkage between marker loci and
  quantitative loci in crosses between inbred lines.
  \emph{Theor. Appl. Genet.} \bold{47}, 35--39. 
}

\author{Karl W Broman, \email{kbroman@jhsph.edu}; Hao Wu,
  \email{hao@jax.org} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
fake.f2 <- calc.genoprob(fake.f2, step=10)
out.2dim <- scantwo(fake.f2, method="hk")
plot(out.2dim)

\dontshow{permo.2dim <- scantwo(fake.f2, method="hk", n.perm=2)
}\dontrun{permo.2dim <- scantwo(fake.f2, method="hk", n.perm=1000)
}apply(permo.2dim,2,quantile,0.95)

# covariates
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=10)
ac <- fake.bc$pheno[,c("sex","age")]
ic <- fake.bc$pheno[,"sex"]
out <- scantwo(fake.bc, method="hk", pheno.col=1,
               addcovar=ac, intcovar=ic)
plot(out)
}

\seealso{ \code{\link[qtl]{plot.scantwo}}, \code{\link[qtl]{summary.scantwo}},
  \code{\link[qtl]{scanone}}, \code{\link[qtl]{max.scantwo}} }

\keyword{models}
