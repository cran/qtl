\name{makeqtl}
\alias{makeqtl}

\title{Make a qtl object}

\description{
  This function takes a cross object and specified chromosome numbers
  and positions and pulls out the genotype probabilities and/or imputed
  genotypes at the nearest pseudomarkers, for later use by the function
  \code{\link[qtl]{fitqtl}}. 
}

\usage{
makeqtl(cross, chr, pos, qtl.name)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Vector indicating the chromosome for each QTL.}
 \item{pos}{Vector (of same length as \code{chr}) indicating the
   positions on the chromosome to be taken. If there's no marker or
   pseudomarker at a position, genotypes for the nearest positions are
   taken.} 
 \item{qtl.name}{The user-specified name for each QTL, used in the
   drop-one-term ANOVA table in \code{\link[qtl]{fitqtl}}.
   If unspecified, the names will be of the form \code{"Chr1@10"} for a
   QTL on Chromsome 1 at 10 cM.} 
}

\details{
  This function will take out the genotype probabilities and imputed
  genotypes if they are present in the input \code{cross} object. If both
  fields are missing in the input object, the function will report an
  error. Before running this function, the user must have first run either
  \code{\link[qtl]{calc.geno}} or \code{\link[qtl]{sim.geno}}.
}

\value{
  An object of class \code{qtl} with the following elements:
  \item{geno}{Imputed genotypes.}
  \item{prob}{Genotype probabilities.}
  \item{chr}{Input vector of chromosome numbers.}
  \item{pos}{Input vector of chromosome positions.}
  \item{n.qtl}{Number of QTLs.}
  \item{n.ind}{Number of individuals.}
  \item{n.gen}{A vector indicating the number of genotypes for each QTL.}
}

\references{
  Sen, S. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Hao Wu, \email{hao@jax.org} }

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp)
}

\seealso{ \code{\link[qtl]{fitqtl}}, \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{sim.geno}} } 

\keyword{models}
