\name{scanone}
\alias{scanone}

\title{Genome scan with single QTL model}

\description{
  Uses interval mapping, ANOVA, Haley-Knott regression, or imputation to
  perform a genome scan, calculating LOD scores under a single QTL model.
  Currently we don't allow covariates, but that will be added
  eventually.  
}

\usage{
scanone(cross, chr, pheno.col=1, method=c("im","anova","hk","imp"),
        start=NULL,maxit=1000, tol=1e-8)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Vector indicating the chromosomes for which LOD scores
   should be calculated.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.}
 \item{method}{Indicates whether to use interval mapping, analysis of
   variance, Haley-Knott regression or imputation.}
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
  \code{\link[qtl]{calc.genoprob}}.  The imputation method uses the
  results of \code{\link[qtl]{sim.geno}}.

  Interval mapping uses the full maximum likelihood method of Lander and
  Botstein (1992).  Analysis of variance (anova) is simply anova using
  the genotypes at marker locations (individuals with missing genotypes
  are discarded).  Haley-Knott regression uses the regression of
  phenotypes on multipoint genotype probabilities.  The imputation
  method uses the pseudomarker algorithm described by Sen and Churchill
  (2001). 

  Individuals with missing phenotypes are dropped.
}

\value{
  A matrix of dimension [n.pos x (n.gen+4)].  The columns are 
  chromosome number, position (in cM; female cM for the 4-way cross),
  LOD score, genotype-specific average phenotypes, and residual SD. 
}

\references{
  CS Haley, SA Knott (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  Heredity 69:315-324.

  ES Lander, D Botstein (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps.  Genetics 121:185-199.

  S Sen, GA Churchill (2001) A statistical framework for quantitative
  trait mapping.  Genetics 159:371-387.

  M Soller, T Brody, A Genizi (1976) On the power of experimental
  designs for the detection of linkage between marker loci and
  quantitative loci in crosses between inbred lines.  Theor Appl Genet
  47:35-39. 
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.bc)
fake.bc$pheno[,2] <- log(fake.bc$pheno[,2])
out.anova <- scanone(fake.bc,pheno.col=2,method="anova")
fake.bc <- calc.genoprob(fake.bc,step=2,off.end=5)
out.im <- scanone(fake.bc,pheno.col=2,method="im")
out.hk <- scanone(fake.bc,pheno.col=2,method="hk")
plot(out.hk,out.anova,lty=1,col=c("red","blue"),chr=c(1,9,11))
plot(out.im,add=TRUE,lty=1,col="black",chr=c(1,9,11))
summary(out.im,3)
}

\seealso{\code{\link[qtl]{plot.scanone}}, \code{\link[qtl]{summary.scanone}},
  \code{\link[qtl]{scanone.perm}}, \code{\link[qtl]{calc.genoprob}} }

\keyword{models}
\keyword{regression}
