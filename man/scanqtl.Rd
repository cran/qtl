\name{scanqtl}
\alias{scanqtl}

\title{General QTL scan}

\description{
  Performs a multiple QTL scan for specified chromosomes and positions
  or intervals, with the possible inclusion of QTL-QTL interactions and/or
  covariates. 
}

\usage{
scanqtl(cross, pheno.col=1, chr, pos, covar=NULL, formula,
        method=c("imp"), incl.markers=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.}
 \item{chr}{Vector indicating the chromosome for each QTL.}
 \item{pos}{List indicating the positions or intervals on the chromosome
   to be scanned.  Each element should be either a single number (for a
   specific position) or a pair of numbers (for an interval).}
 \item{covar}{A data frame of covariates.}
 \item{formula}{An object of class \code{\link[base]{formula}}
   indicating the model to be fitted.  QTLs are indicated as \code{Q1},
   \code{Q2}, etc.  Covariates are indicated by their names in \code{covar}.}
 \item{method}{Indicates whether to use the EM algorithm or  
   imputation. (Only imputation is implemented at this point.)}
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.}
}

\details{
  The formula is used to specified the model to be fitted. In the
  formula, user has to use \code{Q1}, \code{Q2}, etc., or \code{q1},
  \code{q2}, etc., to represent the QTLs, and the column names in the
  covariate data frame to represent the covariates.
    
  Only the interaction terms need to be specifed in the formula. The
  main effects of all input QTLs (as specified by chr and pos) and
  covariates (as specifed by covar) will be included by default. For
  example, if the formula is \code{y~Q1*Q2*Sex}, there are three
  elements in input chr and pos and Sex is one of the column names for
  input cov, the formula used in genome scan will be
  \code{y ~ Q1 + Q2 + Q3 + Sex + Q1:Q2 + Q1:Sex + Q2:Sex + Q1:Q2:Sex}. 

  The input pos is a list or vector to specify the position/range of the
  input chromosomes to be scanned. If it is a vector, it gives the
  specified position on the chromosomes. If it is a list, it will
  contain either the certain positions or a range on chromosomes. For
  example, the input chr is c(1, 6, 13). If pos is c(9.8, 34.0, 18.6),
  it means to scan chromosome 1 at 9.8cM, chromosome 6 at 34cM and 
  chromosome 13 at 18.6cM. If pos is list(c(5,15), c(30,36), 18), it
  means to scan chromosome 1 from 5cM to 15cM, chromosome 6 from 30cM to
  36cM, fix chromosome 13 at 18cM.
}
 
\value{
  An object of class \code{scanqtl}. It is a multi-dimensional
  array of LOD scores, with the number of dimension equal to the number
  of QTLs specifed.} 

\references{
  Sen, S. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Hao Wu, \email{hao@jax.org} }

\examples{
data(fake.f2)

# take out several QTLs
qc <- c(1, 8, 13)
fake.f2 <- subset(fake.f2, chr=qc)

# imputate genotypes
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)
fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)}

# 2-dimensional genome scan with additive 3-QTL model
pos <- list(c(15,35), c(45,65), 28)
result <- scanqtl(fake.f2, pheno.col=1, chr=qc, pos=pos,
                  formula=y~Q1+Q2+Q3)

# image of the results
# chr locations
chr1 <- as.numeric(matrix(unlist(strsplit(colnames(result),"@")),
                   ncol=2,byrow=TRUE)[,2])
chr8 <- as.numeric(matrix(unlist(strsplit(rownames(result),"@")),
                   ncol=2,byrow=TRUE)[,2])
# image plot
image(chr1, chr8, t(result), las=1, col=rev(rainbow(256,start=0,end=2/3)))

# do the same, allowing the QTLs on chr 1 and 13 to interact
result2 <- scanqtl(fake.f2, pheno.col=1, chr=qc, pos=pos,
                   formula=y~Q1+Q2+Q3+Q1:Q3)
# image plot
image(chr1, chr8, t(result2), las=1, col=rev(rainbow(256,start=0,end=2/3)))
}

\seealso{ \code{\link[qtl]{fitqtl}}, \code{\link[qtl]{makeqtl}} }

\keyword{models}
