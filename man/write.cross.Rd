\name{write.cross}
\alias{write.cross}

\title{Write data for a QTL experiment to a file}

\description{
  Data for a QTL experiment is written to a file (or files) in either
  comma-delimited or Mapmaker format.
}

\usage{
write.cross(cross, format=c("csv","mm"), filestem="data", chr, digits=5)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
    \code{\link[qtl]{read.cross}} for details.}
  \item{format}{Specifies whether to write the data in comma-delimited
    or Mapmaker format.}
  \item{filestem}{A character string giving the first part of the output
    file names (the bit before the dot).}
  \item{chr}{A vector specifying for which chromosomes data should be
    written.}
  \item{digits}{Number of digits to which phenotype values should be
    rounded.}  
}

\value{
  Data is written to a single file (comma-delimited format) or two files
  (Mapmaker format).  
}

\details{
  Data is passed to either \code{\link[qtl]{write.cross.csv}} or
  \code{\link[qtl]{write.cross.mm}}, according to the specified
  \code{format} argument.
}

\examples{
\dontrun{data{fake.bc}}
\dontrun{write.cross(fake.bc,"csv","Data/fakebc",c(5,13))}
}

\author{Karl W Broman, 
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{\code{\link[qtl]{read.cross}}, \code{\link[qtl]{write.cross.mm}},
         \code{\link[qtl]{write.cross.csv}} }

\keyword{file}
