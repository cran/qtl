\name{read.cross.csv}
\alias{read.cross.csv}

\title{Read data for a QTL experiment in comma-delimited format}

\description{
  Data for a QTL experiment is read from a file and converted
  into an object of class \code{cross}.  
}

\usage{
read.cross.csv(dir,file,sep=",",na.strings="-",genotypes=c("A","H","B","C","D"))
}			    

\arguments{
  \item{dir}{Character string specifying the directory containing the
    input files.}
  \item{file}{A comma-delimited file containing all of the data, in a
    very rigid format.  The first line should contain the phenotype
    names followed by the marker names.  The next line should contain
    blanks in the phenotype columns, followed by chromosome identifiers
    for each marker in all other columns.  (\bold{Note}: if a chromosome
    has the identifier \code{"X"} or \code{"x"}, it is assumed to be the
    X chromosome.)  The third line should contain blanks in the
    phenotype columns, followed by marker positions, in cM.  Subsequent
    lines should give the phenotype data followed by the genotype data,
    with one line for each individual.}
  \item{sep}{The field separator.  This is generally \code{","}, but
    could be any other character, as long as that character does not
    appear in any of the records.}
  \item{na.strings}{A vector of strings which are to be interpreted as
    missing values.  These are interpreted globally for the entire file,
    so missing value codes in phenotypes must not be valid genotypes,
    and vice versa.}
  \item{genotypes}{A vector of character strings specifying the genotype
    codes used in the file.  Generally this is a vector of length 5,
    with the elements corresponding to AA, AB, BB, not AA (i.e., AB or
    BB), and not BB (ie, AB or BB).  \bold{Note}: Pay careful attention
    to the third and fourth of these; the order of these can be
    confusing!  These five character strings are converted to
    (1,2,3,5,4) in the \code{cross} object that is returned.}
}
\value{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}} for details.
}

%\details{
%}

\examples{
\dontrun{cross <- read.cross("csv",dir="../Data",file="dataset.csv")}
}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{ \code{\link[qtl]{read.cross}},
  \code{\link[qtl]{read.cross.karl}},
  \code{\link[qtl]{read.cross.mm}},
  \code{\link[qtl]{read.cross.csv}} }

\keyword{file}
