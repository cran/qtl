\name{write.cross.csv}
\alias{write.cross.csv}

\title{Write data in comma-delimited format}

\description{
  Data for a QTL experiment is written to a .csv file.
}

\usage{
write.cross.csv(cross, filestem="data", chr, digits=5)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
    \code{\link[qtl]{read.cross}} for details.}
  \item{filestem}{A character string giving the first part of the output
    file names (the bit before the dot).}
  \item{chr}{A vector specifying for which chromosomes data should be
    written.}
  \item{digits}{Number of digits to which phenotype values should be
    rounded.}  
}

\value{
  Data is written to a comma-delimited file, \code{"filestem.csv"},
  where \code{"filestem"} is to be replaced by the corresponding
  argument to the function.  See \code{\link[qtl]{read.cross.csv}} for a
  description of the format.
}

\examples{
\dontrun{data{fake.bc}}
\dontrun{write.cross.csv(fake.bc,"Data/fakebc",c(5,13))}
\dontrun{write.cross(fake.bc,"csv","Data/fakebc",c(5,13))}
}

\author{Karl W Broman, 
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{\code{\link[qtl]{read.cross}}, \code{\link[qtl]{write.cross}},
         \code{\link[qtl]{read.cross.csv}} }

\keyword{file}
