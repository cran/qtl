\name{read.cross.karl}
\alias{read.cross.karl}

\title{Read data for a QTL experiment in Karl Broman's format}

\description{
  Data for a QTL experiment is read from a set of files and converted
  into an object of class \code{cross}.
}

\usage{
read.cross.karl(dir,genfile,mapfile,phefile)
}			    

\arguments{
  \item{dir}{Character string specifying the directory containing the
    input files.}
  \item{genfile}{File containing a (n.ind x n.mar) matrix of
    genotypes, coded 1/2/3/4/5/0 for AA/AB/BB/not BB/not AA/missing,
    one individual per line.}

  \item{phefile}{File containing a (n.ind x n.phe) matrix of
    phenotypes (one individual per line) with "-" denoting
    missing. The first line in the file should give the phenotype names.}

  \item{mapfile}{File containing the map information, in the following
     format \cr
    \code{n.chr}\cr
    \code{n.mar(1) rf(1,1) rf(1,2) \ldots rf(1,n.mar(1)-1)}\cr
    \code{mar.name(1,1)}\cr 
    \code{mar.name(1,2)}\cr 
    \code{\ldots}\cr 
    \code{mar.name(1,n.mar(1))}\cr 
    \code{n.mar(2)}\cr
    \code{\ldots}\cr 
    \code{etc.} }
  }
}
\value{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}} for details.
}

\details{
  We use the marker names to infer the chromosome names (assuming that
  most marker names are like D1M120 or D19M99); if the markers look like
  DXM*, we assume it is the X chromosome.  If all genotypes are <= 2, we
  assume it is a backcross, otherwise it is assumed to be an
  intercross. We weren't thinking about 4-way crosses, and so this
  format doesn't work for that type of cross.
}

\examples{
\dontrun{cross2 <- read.cross("karl",dir="Data", genfile="gen.txt",
		     phefile="phe.txt", mapfile="map.txt")}
}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{ \code{\link[qtl]{read.cross}},
  \code{\link[qtl]{read.cross.mm}}, \code{\link[qtl]{read.cross.gary}},
  \code{\link[qtl]{read.cross.csv}} }

\keyword{file}
