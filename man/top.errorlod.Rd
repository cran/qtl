\name{top.errorlod}
\alias{top.errorlod}

\title{List genotypes with highest error LOD values}

\description{
  Pulls from the output of \code{\link[qtl]{calc.errorlod}} those genotypes
  with error LOD scores above a specified cutoff.
}

\usage{
top.errorlod(cross, chr, cutoff=2, msg=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector specifying for which chromosomes the error LOD
   scores should be inspected.}
 \item{cutoff}{Only those genotypes with error LOD scores above this
   cutoff will be listed.}
 \item{msg}{If TRUE, print a message if there are no apparent errors.}
}

\value{
  A matrix with 4 columns, whose rows correspond to the genotypes that
  are possibly in error.  The four columns give the chromosome number,
  individual number, marker name, and error LOD score.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2)
fake.f2 <- calc.errorlod(fake.f2)
top.errorlod(fake.f2,c(5,13),3)
}

\seealso{\code{\link[qtl]{calc.errorlod}},
  \code{\link[qtl]{plot.errorlod}}, \code{\link[qtl]{find.errors}},
  \code{\link[qtl]{plot.errors}} }

\keyword{manip}
