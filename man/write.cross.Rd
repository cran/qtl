\name{write.cross}
\alias{write.cross}

\title{Write data for a QTL experiment to a file}

\description{
  Data for a QTL experiment is written to a file (or files).
}

\usage{
write.cross(cross, format=c("csv","mm","qtlcart","gary"),
            filestem="data", chr, digits=5)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
    \code{\link[qtl]{read.cross}} for details.}
  \item{format}{Specifies whether to write the data in comma-delimited,
    Mapmaker, QTL Cartographer, or Gary Churchill's format.}
  \item{filestem}{A character string giving the first part of the output
    file names (the bit before the dot).  In Windows, use forward
    slashes (\code{"/"}) or double backslashes (\code{"\\\\"}) to
    specify directory trees.} 
  \item{chr}{A vector specifying for which chromosomes genotype data
    should be written.}
  \item{digits}{Number of digits to which phenotype values should be
    rounded.}
}

%\value{
%  Data is written to a single file (comma-delimited format) or two files
%  (Mapmaker format).
%}

\details{
  Comma-delimited format: a single csv file is created in the format
  read by \code{\link[qtl]{read.cross}}.

  Mapmaker format: Data is written to two files.
  Suppose \code{filestem="file"}.  Then \code{"file.raw"} will contain
  the genotype and phenotype data, and \code{"file.prep"} will contain
  the necessary code for defining the chromosome assignments, marker
  order, and inter-marker distances.

  QTL Cartographer format: Data is written to two files.    Suppose
  \code{filestem="file"}.  Then \code{"file.cro"} will contain 
  the genotype and phenotype data, and \code{"file.map"} will contain
  the genetic map information.  Note that cross types are converted to
  QTL Cartographer cross types as follows: riself to RF1, risib to RF2,
  bc to B1 and f2 to RF2. 

  Gary's format: Data is written to six files. They are:
  \code{"geno.data"} - genotype data; 
  \code{"pheno.data"} - phenotype data; 
  \code{"chrid.dat"} - the chromosome identifier for each marker; 
  \code{"mnames.txt"} - the marker names; 
  \code{"markerpos.txt"} - the marker positions; 
  \code{"pnames.txt"} - the phenotype names
}


\examples{
\dontrun{data(fake.bc)

# comma-delimited format
write.cross(fake.bc, "csv", "Data/fakebc", c(1,5,13))

# Mapmaker format
write.cross(fake.bc, "mm", "Data/fakebc", c(1,5,13))

# QTL Cartographer format
write.cross(fake.bc, "qtlcart", "Data/fakebc", c(1,5,13))

# Gary's format
write.cross(fake.bc, "gary", c(1,5,13))} }

\author{Karl W Broman, \email{kbroman@jhsph.edu}; Hao Wu; Brian S. Yandell}

\seealso{ \code{\link[qtl]{read.cross}} }

\keyword{IO}
