\name{read.cross}
\alias{read.cross}

\title{Read data for a QTL experiment}

\description{
  Data for a QTL experiment is read from a set of files and converted
  into an object of class \code{cross}.  The comma-delimited format
  (\code{csv}) is recommended.  All formats require chromosome
  assignments for the genetic markers, and assume that markers are in
  their correct order.
}

\usage{
read.cross(format=c("csv","mm","gary","karl"), dir=".", file,
           genfile, mapfile, phefile, chridfile, mnamesfile,
           pnamesfile, sep=",", na.strings="-",
           genotypes=c("A","H","B","C","D"), estimate.map=FALSE)
}

\arguments{
  \item{format}{Specifies the format of the data.}
  \item{dir}{Directory in which the data files will be found.  In
    Windows, use forward slashes (\code{"/"}) or double backslashes
    (\code{"\\\\"}) to specify directory trees.}
  \item{file}{The main imput file for formats \code{csv} and \code{mm}.}
  \item{genfile}{File with genotype data (formats \code{karl} and
    \code{gary} only).}
  \item{mapfile}{File with marker position information (all formats
    except \code{csv}).}
  \item{phefile}{File with phenotype data (formats \code{karl} and
    \code{gary} only).}
  \item{chridfile}{File with chromosome ID for each marker (\code{gary}
    format only).} 
  \item{mnamesfile}{File with marker names (\code{gary} format only).}
  \item{pnamesfile}{File with phenotype names (\code{gary} format
    only).}
  \item{sep}{The field separator (\code{csv} format only).  This is
    generally \code{","}, but could be any other character (such as
    \code{"\\t"} for tab), as long as that character does not appear in
    any of the records.}  
  \item{na.strings}{A vector of strings which are to be interpreted as
    missing values (\code{csv} format only).  These are interpreted
    globally for the entire file, so missing value codes in phenotypes
    must not be valid genotypes, and vice versa.}
  \item{genotypes}{A vector of character strings specifying the genotype
    codes (\code{csv} format only).  Generally this is a vector of
    length 5, with the elements corresponding to AA, AB, BB, not AA 
    (i.e., AB or  BB), and not BB (ie, AB or BB).  \bold{Note}: Pay
    careful attention to the third and fourth of these; the order of
    these can be confusing!  These five character strings are converted
    to (1,2,3,5,4) in the \code{cross} object that is returned.} 
  \item{estimate.map}{For formats \code{csv} and \code{mm} only: if TRUE
    and marker positions are not included in the input files, the
    genetic map is estimated using the function
    \code{\link[qtl]{est.map}}.}   
}

\value{
  An object of class \code{cross}, which is a list with two components:

  \item{geno}{This is a list with elements corresponding to
    chromosomes.  \code{names(geno)} contains the names of the
    chromsomes.  Each chromosome is itself a list, and is given class
    \code{A} or \code{X} according to whether it is autosomal
    or the X chromosome. \cr \cr
    
    There are two components for each chromosome: \code{data}, a matrix
    whose rows are individuals and whose columns are markers, and
    \code{map}, either a vector of marker positions (in cM) or a matrix
    of dim (\code{2 x n.mar}) where the rows correspond to marker
    positions in female and male genetic distance, respectively.  \cr \cr

    The genotype data for a backcross is coded as follows: NA = missing,
    1 = AA, 2 = AB.  \cr \cr

    For an F2 intercross, the coding is NA = missing, 1 = AA, 2 = AB, 3
    = BB, 4 = not BB (ie AA or AB; D in mapmaker/qtl), 5 = not AA (ie AB
    or BB; C in mapmaker/qtl). \cr \cr

    For a 4-way cross, the mother and father are assumed to have
    genotypes AB and CD, respectively.  The genotype data for the
    progeny is assumed to be phase-known, with the following coding
    scheme: NA = missing, 1 = AC, 2 = BC, 3 = AD, 4 = BD, 5 = A = AC or AD,
    6 = B = BC or BD, 7 = C = AC or BC, 8 = D = AD or BD, 9 = AC or BD,
    10 = AD or BC.
  }

  \item{pheno}{data.frame of size (\code{n.ind x n.phe}) containing the
    phenotypes.} 
}

\details{
  The available formats are comma-delimited (\code{csv}), Mapmaker
  (\code{mm}), Gary Churchill's format (\code{gary}) and Karl Broman's
  format (\code{karl}).  The required files and their specification for
  each format appears below.  The comma-delimited format is recommended.
  Note that these formats work only for backcross and intercross data.
}


\section{CSV format}{
  The input file is a text file with a specified field delimiter
  (\code{sep}) (a comma is recommended).

  The first line should contain the phenotype names followed by the
  marker names.  \bold{At least one phenotype must be included}; for
  example, include a numerical index for each individual.

  The second line should contain blanks in the phenotype columns,
  followed by chromosome identifiers for each marker in all other
  columns. If a chromosome has the identifier \code{X} or \code{x}, it
  is assumed to be the X chromosome; otherwise, it is assumed to be an
  autosome.

  An optional third line should contain blanks in the phenotype
  columns, followed by marker positions, in cM. 

  Subsequent lines should give the data, with one line for each
  individual, and with phenotypes followed by genotypes.  If possible,
  phenotypes are made numeric; otherwise they are converted to factors.

  The cross is determined to be a backcross if only the first two elements
  of the \code{genotypes} string are found; otherwise, it is assumed to
  be an intercross.
}

\section{Mapmaker format}{
  This format requires two files.  The so-called rawfile, specified by
  the argument \code{file}, contains the genotype and phenotype
  data. Rows beginning with the symbol \code{#} are ignored.  The first
  line should be either \code{data type f2 intercross} or
  \code{data type f2 backcross}.  The second line should begin with
  three numbers indicating the numbers of individuals, markers and
  phenotypes in the file.  This line may include the word \code{symbols}
  followed by symbol assignments (see the documentation for mapmaker,
  and cross your fingers).  The rest of the lines give genotype data
  followed by phenotype data, with marker and phenotype names always
  beginning with the symbol \code{*}.
    
  A second file, giving marker positions and chromosome assignments, is
  not standard for the Mapmaker software, but is necessary for us; it is
  specified with the argument \code{mapfile}.  (For this format, if
  \code{genfile} is specified by not \code{mapfile}, we assume that
  \code{genfile} is the file to use.)  The file should contain two or
  three columns separated by white space and with no header row.
  The first column gives the chromosome assignments.  The second column
  gives the marker names, with markers listed in the order along the
  chromosomes.  An optional third column lists the map positions of the
  markers.

  If a chromosome has the identifier \code{X} or \code{x}, it is assumed
  to be the X chromosome; otherwise, it is assumed to be an autosome.
}

\section{Gary format}{
  This format requires the six files; all files have default names, and
  so need not be specified if the default name is used. 

  \code{genfile} (default = \code{"geno.dat"}) contains the genotype
  data.  The file contains one line per individual, with genotypes for
  the set of markers separated by white space.  Missing values are
  coded as 9, and genotypes are coded as 0/1/2 for AA/AB/BB.

  \code{mapfile} (default = \code{"markerpos.txt"}) contains the two
  columns with no header row: the marker names in the first columns and
  their cM position in the second column.
  
  \code{phefile} (default = \code{"pheno.dat"}) contains the phenotype
  data, with one row for each mouse and one column for each phenotype.
  There should be no header row, and missing values are coded as
  \code{"-"}. 

  \code{chridfile} (default = \code{"chrid.dat"}) contains the
  chromosome identifier for each marker.  

  \code{mnamesfile} (default = \code{"mnames.txt"}) contains the marker
  names.

  \code{pnamesfile} (default = \code{"pnames.txt"}) contains the names
  of the phenotypes.
}


\section{Karl format}{
  This format requires three files; all files have default names, and so
  need not be specified if the default name is used. 

  \code{genfile} (default = \code{"gen.txt"}) contains the genotype
  data.  The file contains one line per individual, with genotypes
  separated by white space.  Missing values are coded 0; genotypes are
  coded as 1/2/3/4/5 for AA/AB/BB/not BB/not AA.
    
  \code{mapfile} (default = \code{"map.txt"}) contains the map
  information, in the following complicated format: \cr \cr
    \code{n.chr} \cr
    \code{n.mar(1) rf(1,1) rf(1,2) \ldots rf(1,n.mar(1)-1)}\cr
    \code{mar.name(1,1)}\cr 
    \code{mar.name(1,2)}\cr 
    \code{\ldots}\cr 
    \code{mar.name(1,n.mar(1))}\cr 
    \code{n.mar(2)}\cr
    \code{\ldots}\cr 
    \code{etc.} \cr

  \code{phefile} (default = \code{"phe.txt"}) contains a matrix of
  phenotypes, with one individual per line.  The first line in the
  file should give the phenotype names.
}

\examples{
\dontrun{# comma-delimited format
dat1 <- read.cross("csv", dir="Mydata", file="mydata.csv")

# Mapmaker format
dat2 <- read.cross("mm", dir="Mydata", file="mydata.raw",
                   mapfile="mydata.map")

# Gary format
dat3 <- read.cross("gary", dir="Mydata", genfile="geno.dat",
                   mapfile="markerpos.txt", phefile="pheno.dat",
                   chridfile="chrid.dat", mnamesfile="mnames.txt",
                   pnamesfile="pnames.txt")

# Karl format
dat4 <- read.cross("karl", dir="Mydata", genfile="gen.txt",
                   phefile="phe.txt", mapfile="map.txt")}
} 

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\seealso{ \code{\link[qtl]{write.cross}}, \code{\link[qtl]{sim.cross}};
  the \code{sampledata} directory in the package distribution contains
  sample data files in all formats except Gary's.
}

\keyword{IO}
