\name{read.cross.mm}
\alias{read.cross.mm}

\title{Read data for a QTL experiment in mapmaker format}

\description{
  Data for a QTL experiment is read from a set of files and converted
  into an object of class \code{cross}.
}

\usage{
read.cross.mm(dir,rawfile,mapfile,estimate.map)
}			    

\arguments{
  \item{dir}{Character string specifying the directory containing the
    input files.}
  \item{rawfile}{Mapmaker rawfile containing the genotype and phenotype
    data. Rows beginning with the symbol \code{#} are ignored.  The
    first line should be either \code{data type f2 intercross} or
    \code{data type f2 backcross}.  The second line should begin with
    three numbers indicating the numbers of individuals, markers and
    phenotypes in the file.  This line may include the word
    \code{symbols} followed by symbol assignments (see the documentation
    for mapmaker, and cross your fingers).  The rest of the lines give
    genotype data followed by phenotype data, with marker and phenotype
    names always beginning with the symbol \code{*}.
  }

  \item{mapfile}{File containing two or three columns separated by white
    space, with \emph{no} header row.  The first column gives the
    chromosome assignments.  The second column gives the marker names,
    with markers listed in the order along the chromosomes.  An optional
    third column lists the map positions of the markers.}

  \item{estimate.map}{If \code{mapfile} does not contain the markers' map
    positions and \code{estimate.map} is \code{TRUE}, the genetic map
    for the cross is estimated using the function \code{\link[qtl]{est.map}}.}

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
\dontrun{cross1 <- read.cross(format="mm", dir="Data",
rawfile="sample.raw", mapfile="sample.map")}		     
}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{ \code{\link[qtl]{read.cross}},
  \code{\link[qtl]{read.cross.karl}},
  \code{\link[qtl]{read.cross.gary}},
  \code{\link[qtl]{read.cross.csv}} }
}

\keyword{file}
