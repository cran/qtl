\name{read.cross.gary}
\alias{read.cross.gary}

\title{Read data for a QTL experiment in Gary Churchill's format}

\description{
  Data for a QTL experiment is read from a set of files and converted
  into an object of class \code{cross}.  
}

\usage{
read.cross.gary(dir,genfile,mnamesfile,chridfile,phefile,pnamesfile,mapfile)
}			    

\arguments{
  \item{dir}{Character string specifying the directory containing the
    input files.}
  \item{genfile}{File containing genotype data.}
  \item{mnamesfile}{File containing a list of marker names.}
  \item{chridfile}{File containing a list of chromosome numbers, one for
    each marker.}
  \item{phefile}{File containing phenotype data.}
  \item{pnamesfile}{File containing a list of phenotype names.}
  \item{mapfile}{File containing map locations for the markers.  There
    may be extraneous records in this file.}
}
\value{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}} for details.
}

%\details{
%}

\examples{
\dontrun{cross <- read.cross("gary",dir="../Data",genfile="geno.dat",
      mnamesfile="mnames.txt",chridfile="chrid.dat",phefile="pheno.dat",
      pnamesfile="pnames.txt",mapfile="markerpos.txt")}
}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{ \code{\link[qtl]{read.cross}},
  \code{\link[qtl]{read.cross.karl}},
  \code{\link[qtl]{read.cross.mm}},
  \code{\link[qtl]{read.cross.csv}} }

\keyword{file}
