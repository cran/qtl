\name{scanone.perm}
\alias{scanone.perm}

\title{Permutation test for genome scan with single QTL model}

\description{
  Performs a permutation test of interval mapping, ANOVA or Haley-Knott
  regression in a genome scan under a single QTL model.
  Currently we don't allow covariates, but that will be added
  eventually.  
}

\usage{
scanone.perm(cross, chr, pheno.col=1, method=c("im","anova","hk"),
             start=NULL,n.perm=1000, maxit=1000, tol=1e-8)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Vector indicating the chromosomes for which LOD scores
   should be calculated.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.}
 \item{method}{Indicates whether to use interval mapping, analysis of
   variance or Haley-Knott regression.}
 \item{n.perm}{Number of permutations to be performed.}
 \item{start}{If \code{NULL}, use the usual starting values; if length
   1, use random initial weights for EM; otherwise, this should be a
   vector of length N+1 (where N is the number of possible genotypes for
   the cross), giving the initial values for EM.}
 \item{maxit}{Maximum number of iterations in the EM algorithm; used
   only in interval mapping.}
 \item{tol}{Tolerance value for determining convergence in the EM
   algorithm; used only in interval mapping.}
}

\details{
  The interval mapping and Haley-Knott regression methods require that
  multipoint genotype probabilities are first calculated using
  \code{\link[qtl]{calc.genoprob}}.  The analysis of variance method
  uses reconstructed genotypes calculated using
  \code{\link[qtl]{argmax.geno}}. 

  Interval mapping uses the full maximum likelihood method of Lander and
  Botstein (1992).  Analysis of variance (anova) is simply anova using
  the genotypes at marker locations reconstructed by the Viterbi
  algorithm.  Haley-Knott regression uses the regression of phenotypes
  on multipoint genotype probabilities.

  Individuals with missing phenotypes are dropped.
}

\value{
  A matrix of dimension [n.perm x (n.chr+1)].  The columns are 
  the maximum LOD score across the whole genome and across each
  respective chromosome, for each of the permutation replicates.
}

\references{
  GA Churchill, RW Doerge (1994) Empirical threshold values for
  quantitative trait mapping.  Genetics 138:963-971.

  CS Haley, SA Knott (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  Heredity 69:315-324.

  ES Lander, D Botstein (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps.  Genetics 121:185-199.

  M Soller, T Brody, A Genizi (1976) On the power of experimental
  designs for the detection of linkage between marker loci and
  quantitative loci in crosses between inbred lines.  Theor Appl Genet
  47:35-39. 
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2,step=1,off.end=5)
output1 <- scanone(fake.f2,method="im")
\dontrun{permo1 <- scanone.perm(fake.f2,method="im",n.perm=1000)}
\testonly{permo1 <- scanone.perm(fake.f2,method="im",n.perm=3)}

data(fake.bc)
fake.bc <- argmax.geno(fake.bc)
output2 <- scanone(fake.bc,chr=c(1,9),method="anova")
\dontrun{permo2 <- scanone.perm(fake.bc,chr=c(1,9),method="anova",n.perm=1000)}
\testonly{permo2 <- scanone.perm(fake.bc,chr=c(1,9),method="anova",n.perm=3)}
threshold <- quantile(permo2[,1],0.95)
}

\seealso{\code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}} }

\keyword{models}
\keyword{regression}
