\name{vbscan}
\alias{vbscan}

\title{Genome scan for trait with some undefined values}

\description{
  Performs a genome scan by interval mapping for a quantitative
  phenotype for which some individuals' phenotype is undefined (for
  example, the size of a lesion, where some individuals have no
  lesion).  The undefined value must correspond to either the largest
  phenotype or the smallest phenotype.  
}

\usage{
vbscan(cross, chr, pheno.col=1, upper=FALSE,
       maxit=1000, tol=1e-8)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Vector indicating the chromosomes for which LOD scores
   should be calculated.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.}
 \item{upper}{If true, the "undefined" phenotype is the maximum observed
   phenotype; otherwise, it is the smallest observed phenotype.}
 \item{maxit}{Maximum number of iterations in the EM algorithm; used
   only in interval mapping.}
 \item{tol}{Tolerance value for determining convergence in the EM
   algorithm; used only in interval mapping.}
}

\details{
  The function requires that multipoint genotype probabilities are first
  calculated using \code{\link[qtl]{calc.genoprob}}.  

  Describe the model here.

  Individuals with missing phenotypes are dropped.
}

\value{
  A matrix of dimension [n.pos x (2*n.gen+6)].  The columns are 
  chromosome number, position (in cM; female cM for the 4-way cross),
  LOD(p,mu), LOD(p), LOD(mu), estimates of the ps, estimates of the mus, 
  and residual SD.
}

\references{
  Boyartchuk VL, Broman KW, Mosher RE, D'Orazio SEF, Starnbach MN,
  Dietrich WF (2001) Multigenic control of Listeria monocytogenes
  susceptibility in mice. Nature Genetics 27:259-260 \cr
  
  Broman KW, Boyartchuk VL, Dietrich WF (2000) Mapping time-to-death
  quantitative trait loci in a mouse cross with high survival
  rates. Technical Report MS00-04, Department of Biostatistics, Johns
  Hopkins University, Baltimore, MD
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(listeria)
listeria$pheno <- log10(listeria$pheno)
listeria <- calc.genoprob(listeria,step=1,off.end=0,err=0)
output <- vbscan(listeria,upper=TRUE)
plot(output)
}

\seealso{\code{\link[qtl]{plot.scanone}}, \code{\link[qtl]{summary.scanone}},
  \code{\link[qtl]{vbscan.perm}}, \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{scanone}} }

\keyword{models}
\keyword{regression}
