\name{read.cross}
\alias{read.cross}

\title{Read data for a QTL experiment}

\description{
  Data for a QTL experiment is read from a set of files and converted
  into an object of class \code{cross}.
}

\usage{
read.cross(format=c("csv","mm","gary","karl"), \dots)
}

\arguments{
  \item{format}{Specifies the format of the data (\code{"csv"},
    \code{"mm"}, \code{"gary"}, or \code{"karl"}).}
    
  \item{\dots}{All other arguments are passed to functions specific to the
    data format.}
}

\value{
  An object of class \code{cross}, which is a list with two components:

  \item{geno}{This is a list with elements corresponding to
    chromosomes.  \cr \cr
    \code{names(geno)} contains the names of the
    chromsomes.  \cr \cr
    Each chromosome is itself a list, and is given class
    \code{A} or \code{X} according to whether it is autosomal
    or the X chromosome. \cr \cr
    There are two components for each chromosome: \cr
    \code{data}, a matrix whose rows are individuals and whose
    columns are markers \cr
    \code{map}, either a vector of marker
    positions (in cM) or a matrix of dim (\code{2 x n.mar}) where the rows
    correspond to marker positions in female and male genetic distance,
    respectively.  \cr \cr
    The genotype data for a backcross is coded as follows: 0 = missing,
    1 = AA, 2 = AB.  \cr \cr
    For an F2 intercross, the coding is 0 = missing, 1 = AA, 2 = AB, 3
    = BB, 4 = not BB (ie AA or AB; D in mapmaker/qtl), 5 = not AA (ie AB
    or BB; C in mapmaker/qtl). \cr \cr
    For a 4-way cross, the mother and father are assumed to have
    genotypes AB and CD, respectively.  The genotype data for the
    progeny is assumed to be phase-known, with the following coding
    scheme: 0 = missing, 1 = AC, 2 = BC, 3 = AD, 4 = BD, 5 = A = AC or AD,
    6 = B = BC or BD, 7 = C = AC or BC, 8 = D = AD or BD, 9 = AC or BD,
    10 = AD or BC. \cr
  }

  \item{pheno}{Matrix of size (\code{n.ind x n.phe}) containing the
    phenotypes} 
}

\examples{
\dontrun{cross1 <- read.cross("karl",dir="Data", genfile="gen.txt",
		     phefile="phe.txt", mapfile="map.txt")}
}

\details{
  This function calls either \code{\link[qtl]{read.cross.csv}}, 
  \code{\link[qtl]{read.cross.mm}}, \code{\link[qtl]{read.cross.gary}},
  or \code{\link[qtl]{read.cross.karl}}, according to the specified 
  format.
}

\author{Karl W Broman, 
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{\code{\link[qtl]{read.cross.csv}}, \code{\link[qtl]{read.cross.gary}}, 
  \code{\link[qtl]{read.cross.mm}}, \code{\link[qtl]{read.cross.karl}},
  \code{\link[qtl]{fake.bc}}, \code{\link[qtl]{fake.f2}},
  \code{\link[qtl]{fake.4way}}, \code{\link[qtl]{listeria}},
  \code{\link[qtl]{hyper}} }

\keyword{file}
