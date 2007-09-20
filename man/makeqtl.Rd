\name{makeqtl}
\alias{makeqtl}

\title{Make a qtl object}

\description{
  This function takes a cross object and specified chromosome numbers
  and positions and pulls out the genotype probabilities or imputed
  genotypes at the nearest pseudomarkers, for later use by the function
  \code{\link[qtl]{fitqtl}}. 
}

\usage{
makeqtl(cross, chr, pos, qtl.name, what=c("draws","prob"))
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
 \item{what}{Indicates whether to pull out the imputed genotypes or the
   genotype probabilities.}
}

\details{
  This function will take out the genotype probabilities and imputed
  genotypes if they are present in the input \code{cross} object. If both
  fields are missing in the input object, the function will report an
  error. Before running this function, the user must have first run either
  \code{\link[qtl]{sim.geno}} (for \code{what="draws"}) or
  \code{\link[qtl]{calc.genoprob}} (for \code{what="prob"}).
}

\value{
  An object of class \code{qtl} with the following elements (though only
  one of \code{geno} and \code{prob} will be included, according to
  whether \code{what} is given as \code{"draws"} or \code{"prob"}):
  \item{\code{geno}}{Imputed genotypes.}
  \item{\code{prob}}{Genotype probabilities.}
  \item{\code{chr}}{Input vector of chromosome numbers.}
  \item{\code{pos}}{Input vector of chromosome positions.}
  \item{\code{n.qtl}}{Number of QTLs.}
  \item{\code{n.ind}}{Number of individuals.}
  \item{\code{n.gen}}{A vector indicating the number of genotypes for each QTL.}
}

\author{Hao Wu}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="draws")
}

\seealso{ \code{\link[qtl]{fitqtl}}, \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{sim.geno}} } 

\keyword{utilities}
