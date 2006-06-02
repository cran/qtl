\name{scanone}
\alias{scanone}

\title{Genome scan with a single QTL model}

\description{
  Genome scan with a single QTL model, with possible allowance for
  covariates, using any of several possible models for the phenotype and
  any of several possible numerical methods.
}

\usage{
scanone(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
        method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
        addcovar=NULL, intcovar=NULL, weights=NULL,
        upper=FALSE, ties.random=FALSE, start=NULL, maxit=4000,
        tol=1e-4, n.perm, n.node=1, verbose)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Vector indicating the chromosomes for which LOD scores
   should be calculated.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This can be a vector of integers; for methods
   \code{"hk"} and \code{"imp"} this can be considerably faster than doing
   them one at a time.}
 \item{model}{The phenotypic model: the usual normal model, a model for
   binary traits, a two-part model or non-parametric.}
 \item{method}{Indicates whether to use the EM algorithm, 
   imputation, Haley-Knott regression, or marker regression.  Not all
   methods are available for all models.  Marker regression is performed
   either by dropping individuals with missing genotypes (\code{"mr"}),
   or by first filling in missing data using a single imputation
   (\code{"mr-imp"}) or by the Viterbi algorithm (\code{"mr-argmax"}).}
 \item{addcovar}{Additive covariates; allowed only for the normal and
   binary models.}
 \item{intcovar}{Interactive covariates (interact with QTL genotype);
   allowed only for the normal and binary models.}
 \item{weights}{Optional weights of individuals.  Should be either NULL
   or a vector of length n.ind containing positive weights.  Used only
   in the case \code{model="normal"}}
 \item{upper}{Used only for the two-part model; if true, the
   ``undefined'' phenotype is the maximum observed phenotype; otherwise,
   it is the smallest observed phenotype.} 
 \item{ties.random}{Used only for the non-parametric ``model;'' if TRUE,
   ties in the phenotypes are ranked at random.  If FALSE, average ranks
   are used and a corrected LOD score is calculated.}
 \item{start}{Used only for the EM algorithm with the normal model and
   no covariates.  If \code{NULL}, use the usual starting values; if
   length 1, use random initial weights for EM; otherwise, this should
   be a vector of length n+1 (where n is the number of possible
   genotypes for the cross), giving the initial values for EM.}
 \item{maxit}{Maximum number of iterations in the EM algorithm; used
   only in interval mapping.}
 \item{tol}{Tolerance value for determining convergence in the EM
   algorithm; used only in interval mapping.}
 \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
 \item{n.node}{Number of nodes to use for parallelization of
   permutations; only if package \code{snow} is available.}
 \item{verbose}{In the case \code{n.perm} is specified, display
   information about the progress of the permutation tests.}
}

\section{Models}{
  \bold{The normal model} is the standard model for QTL mapping.  The
  residual phenotypic variation is assumed to follow a normal
  distribution, and analysis is analogous to linear regression. 

  \bold{The binary model} is for the case of a binary phenotype, which
  must have values 0 and 1.  The proportions of 1's in the different
  genotype groups are compared.  Currently only methods \code{em} and
  \code{mr} are available for this model.

  \bold{The two-part model} is appropriate for the case of a spike in
  the phenotype distribution (for example, metastatic density when many
  individuals show no metastasis, or survival time following an
  infection when individuals may recover from the infection and fail to
  die).  The two-part model was described by Broman et al. (2000) and
  Boyartchuk et al. (2001).  Individuals with QTL genotype \eqn{g} have
  probability \eqn{p_g}{p[g]} of having an undefined phenotype (the
  spike), while if their phenotype is defined, it comes from a normal
  distribution with mean \eqn{\mu_g}{mu[g]} and common standard
  deviation \eqn{\sigma}{s}. Three LOD scores are calculated:
  LOD(\eqn{p,\mu}{p,mu}) is for the test of the hypothesis that \eqn{p_g
    = p}{p[g] = p} and \eqn{\mu_g = \mu}{mu[g] = mu}.  LOD(\eqn{p}) is
  for the test that \eqn{p_g = p}{p[g] = p} while the \eqn{\mu_g}{mu[g]}
  may vary. LOD(\eqn{\mu}{mu}) is for the test that \eqn{\mu_g =
    \mu}{mu[g] = mu} while the \eqn{p_g}{p[g]} may vary. 

  \bold{With the non-parametric ``model''}, an extension of the
  Kruskal-Wallis test is used; this is similar to the method described
  by Kruglyak and Lander (1995).  In the case of incomplete genotype
  information (such as at locations between genetic markers), the
  Kruskal-Wallis statistic is modified so that the rank for each
  individual is weighted by the genotype probabilities, analgous to
  Haley-Knott regression.  For this method, if the argument
  \code{ties.random} is TRUE, ties in the phenotypes are assigned random
  ranks; if it is FALSE, average ranks are used and a corrected LOD
  score is calculate.  Currently the \code{method} argument is ignored
  for this model.
}

\section{Methods}{
  \bold{\code{em}}: maximum likelihood is performed via the
  EM algorithm (Dempster et al. 1977), first used in this context by
  Lander and Botstein (1989).

  \bold{\code{imp}}: multiple imputation is used, as described by Sen
  and Churchill (2001).

  \bold{\code{hk}}: Haley-Knott regression is used (regression of the
  phenotypes on the multipoint QTL genotype probabilities), as described
  by Haley and Knott (1992).

  \bold{\code{mr}}: Marker regression is used.  Analysis is performed
  only at the genetic markers, and individuals with missing genotypes
  are discarded.
}

\section{Covariates}{
Covariates are allowed only for the normal model, in which case the
model is \eqn{y = \beta_q + A \gamma + Z \delta_q + \epsilon}{y = b[q] +
  A g + Z d[q] + e} where \emph{q} is the unknown QTL genotype, \emph{A}
is a matrix of additive covariates, and \emph{Z} is a matrix of
covariates that interact with the QTL genotype.  The columns of \emph{z}
are forced to be contained in the matrix \emph{A}.

The LOD score is calculated comparing the likelihood of the above
model to that of the null model \eqn{y = \mu + A \gamma + \epsilon}{y =
  m + A g + e}.

Covariates must be numeric matrices.  Individuals with any missing
covariates are discarded.  
}

\section{X chromosome}{
The X chromosome must be treated specially in QTL mapping.

If both males and females are included, male hemizygotes are allowed
to be different from female homozygotes.  Thus, in a backcross, we will
fit separate means for the genotype classes AA, AB, AY, and BY.  In such
cases, sex differences in the phenotype could cause spurious linkage to
the X chromosome, and so the null hypothesis must be changed to allow
for a sex difference in the phenotype.


\tabular{lllccc}{
\bold{BC}
 \tab \tab \bold{Sexes} \tab \bold{Null} \tab \bold{Alternative} \tab \bold{df} \cr
 \tab \tab both sexes   \tab sex         \tab AA/AB/AY/BY        \tab 2         \cr 
 \tab \tab all female   \tab grand mean  \tab AA/AB              \tab 1         \cr
 \tab \tab all male     \tab grand mean  \tab AY/BY              \tab 1         \cr
 \tab \tab              \tab             \tab                    \tab           \cr
 
\bold{F2}
 \tab \bold{Direction} \tab \bold{Sexes}\tab \bold{Null}\tab \bold{Alternative} \tab \bold{df} \cr
 \tab Both             \tab both sexes  \tab femaleF/femaleR/male \tab AA/ABf/ABr/BB/AY/BY \tab 3         \cr
 \tab                  \tab all female  \tab pgm        \tab AA/ABf/ABr/BB      \tab 2         \cr
 \tab                  \tab all male    \tab grand mean \tab AY/BY              \tab 1         \cr
 \tab Forward          \tab both sexes  \tab sex        \tab AA/AB/AY/BY        \tab 2         \cr
 \tab                  \tab all female  \tab grand mean \tab AA/AB              \tab 1         \cr
 \tab                  \tab all male    \tab grand mean \tab AY/BY              \tab 1         \cr
 \tab Backward         \tab both sexes  \tab sex        \tab AB/BB/AY/BY        \tab 2         \cr
 \tab                  \tab all female  \tab grand mean \tab AB/BB              \tab 1         \cr
 \tab                  \tab all male    \tab grand mean \tab AY/BY              \tab 1         \cr
}


}

\details{
  Use of the EM algorithm or Haley-Knott regression require that
  multipoint genotype probabilities are first calculated using
  \code{\link[qtl]{calc.genoprob}}.  The imputation method uses the
  results of \code{\link[qtl]{sim.geno}}.

  Individuals with missing phenotypes are dropped.

  In the case that \code{n.perm} is not missing, so that a permutation
  test is performed, the R function \code{scanone} is called repeatedly.
 
  See further details on the models, the methods and the use of
  covariates below.
}

\value{
  If \code{n.perm} is missing, the function returns a data.frame whose
  first two columns contain the chromosome IDs and cM positions.  The
  third column contains the LOD score.  In the case of the two-part
  model, the third column is LOD(\eqn{p,\mu}{p,mu}), while the fourth
  and fifth columns are LOD(\eqn{p}) and LOD(\eqn{\mu}{mu}).  The data
  frame is given class \code{"scanone"} and
  attributes  \code{"model"}, \code{"method"} and \code{"type"} (the
  latter is the type of cross analyzed). \cr 

  If \code{n.perm} is specified, the function returns either a vector of
  length \code{n.perm}, containing the maximum LOD scores, genome-wide,
  for the permutation replicates. In the case of the two-part model, the
  return value is a matrix of size \code{n.perm x 3}, with columns
  corresponding to the three different LOD scores.
}

\references{
  Boyartchuk V. L., Broman, K. W., Mosher, R. E., D'Orazio
  S. E. F., Starnbach, M. N. and Dietrich, W. F. (2001) Multigenic
  control of \emph{Listeria monocytogenes} susceptibility in
  mice. \emph{Nature Genetics} \bold{27}, 259--260.
  
  Broman, K. W., Boyartchuk, V. L. and Dietrich, W. F. (2000) Mapping
  time-to-death quantitative trait loci in a mouse cross with high
  survival rates. Technical Report MS00-04, Department of Biostatistics,
  Johns Hopkins University, Baltimore, MD.

  Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for
  quantitative trait mapping.  \emph{Genetics} \bold{138}, 963--971.

  Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum
  likelihood from incomplete data via the EM algorithm.  \emph{J. Roy.
    Statist. Soc.} B, \bold{39}, 1--38.

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Kruglyak, L. and Lander, E. S. (1995) A nonparametric approach for
  mapping quantitative trait loci.  \emph{Genetics} \bold{139},
  1421--1428. 

  Lander, E. S. and Botstein, D. (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps.  \emph{Genetics}
  \bold{121}, 185--199.

  Sen, S. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.

  Soller, M., Brody, T. and Genizi, A. (1976) On the power of experimental
  designs for the detection of linkage between marker loci and
  quantitative loci in crosses between inbred lines.
  \emph{Theor. Appl. Genet.} \bold{47}, 35--39. 

  Xu, S., and Atchley, W.R. (1996) Mapping quantitative trait loci for
  complex binary diseases using line crosses. \emph{Genetics}
  \bold{143}, 1417--1424.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu}; Hao Wu,
  \email{hao@jax.org} }

\examples{
###################
# Normal Model
###################
data(hyper)
\dontshow{hyper <- subset(hyper,chr=1:4,ind=1:100)}
# Genotype probabilities for EM and H-K
hyper <- calc.genoprob(hyper, step=2.5)
out.em <- scanone(hyper, method="em")
out.hk <- scanone(hyper, method="hk")

# Summarize results: peaks above 3
summary(out.em, thr=3)
summary(out.hk, thr=3)

# Plot the results
plot(out.hk, out.em)
plot(out.hk, out.em, chr=c(1,4), lty=1, col=c("blue","black"))

# Imputation; first need to run sim.geno
# Do just chromosomes 1 and 4, to save time
hyper.c1n4 <- sim.geno(subset(hyper, chr=c(1,4)),
                       step=2.5, n.draws=8)
out.imp <- scanone(hyper.c1n4, method="imp")
summary(out.imp, thr=3)

# Plot all three results
plot(out.imp, out.hk, out.em, chr=c(1,4), lty=1,
     col=c("red","blue","black"))

# Permutation tests
\dontrun{permo <- scanone(hyper, method="hk", n.perm=1000)
}\dontshow{permo <- scanone(hyper, method="hk", n.perm=3)
}quantile(permo, 0.95)


###################
# Non-parametric
###################
out.np <- scanone(hyper, model="np")
summary(out.np, thr=3)

# Plot with previous results
plot(out.np, chr=c(1,4), lty=1, col="green")
plot(out.imp, out.hk, out.em, chr=c(1,4), lty=1,
     col=c("red","blue","black"), add=TRUE)

###################
# Two-part Model
###################
data(listeria)
\dontshow{listeria <- subset(listeria,chr=c(1,5,13))}
listeria <- calc.genoprob(listeria,step=2.5)
out.2p <- scanone(listeria, model="2part", upper=TRUE)
summary(out.2p, thr=5)

# Plot all three LOD scores together
plot(out.2p, out.2p, out.2p, lodcolumn=c(2,3,1), lty=1, chr=c(1,5,13),
     col=c("red","blue","black"))

# Permutation test
\dontrun{permo <- scanone(listeria, model="2part", upper=TRUE,
                 n.perm=1000)
}\dontshow{permo <- scanone(listeria, model="2part", upper=TRUE,
                 n.perm=3)
}apply(permo, 2, quantile, 0.95)

###################
# Binary model
###################
listeria <- subset(listeria, ind=!is.na(listeria$pheno[,1]))
listeria$pheno[,2] <- rep(0,nind(listeria))
listeria$pheno[listeria$pheno[,1]==264,2] <- 1
out.bin <- scanone(listeria, pheno.col=2, model="binary")
summary(out.bin, thr=3)

# Plot LOD for binary model with LOD(p) from 2-part model
plot(out.bin, out.2p, lodcolumn=c(1,2), lty=1, col=c("black", "red"),
     chr=c(1,5,13))

# Permutation test
\dontrun{permo <- scanone(listeria, pheno.col=2, model="binary",
                 n.perm=1000)
}\dontshow{permo <- scanone(listeria, pheno.col=2, model="binary",
                 n.perm=3)
}quantile(permo, 0.95)

###################
# Covariates
###################
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=c(2,5,10))}
plot(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=2.5)
# genome scans without covariates
out.nocovar <- scanone(fake.bc)
# genome scans with covariates
ac <- fake.bc$pheno[,c("sex","age")]
ic <- fake.bc$pheno[,"sex"]
out.covar <- scanone(fake.bc, pheno.col=1,
                     addcovar=ac, intcovar=ic)
summary(out.nocovar, thr=3)
summary(out.covar, thr=3)
plot(out.covar, out.nocovar, chr=c(2,5,10))
}

\seealso{ \code{\link[qtl]{plot.scanone}}, 
  \code{\link[qtl]{summary.scanone}}, \code{\link[qtl]{scantwo}},
  \code{\link[qtl]{calc.genoprob}}, \code{\link[qtl]{sim.geno}},
  \code{\link[qtl]{max.scanone}} } 

\keyword{models}
