\name{write.cross.mm}
\alias{write.cross.mm}

\title{Write data in mapmaker format}

\description{
  Data for a QTL experiment is written to .raw and .prep files, suitable
  for importing into the mapmaker/exp and mapmaker/qtl programs.
}

\usage{
write.cross.mm(cross, filestem="data", chr, digits=5)
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
  Data is written to two files, in mapmaker format.  Suppose
  \code{filestem="file"}.  Then "file.raw" will contain the genotype and
  phenotype data, and "file.prep" will contain the necessary code for
  defining the chromosome assignments, marker order, and inter-marker
  distances. 
}

\examples{
\dontrun{data{fake.bc}}
\dontrun{write.cross.mm(fake.bc,"Data/fakebc",c(5,13))}
}

\author{Karl W Broman, 
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{\code{\link[qtl]{read.cross}}, \code{\link[qtl]{read.cross.mm}} }

\keyword{file}
