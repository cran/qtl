\name{plot.errors}
\alias{plot.errors}

\title{Plot grid of likely genotyping errors}

\description{
  Plots a grid displaying the output of the \code{\link[qtl]{find.errors}}
  function.
}

\usage{
plot.errors(x, chr, \dots)
}
\arguments{
 \item{x}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector specifying for which chromosomes the calculations
   should be performed.}
 \item{\dots}{Ignored at this point.}
}

\details{
  The most likely sequence of underlying genotypes is calculated in
  \code{\link[qtl]{argmax.geno}}.  The function
  \code{\link[qtl]{plot.missing}} is used to make the plot.
}

\value{
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
output <- find.errors(fake.f2)
plot.errors(fake.f2)
}

\seealso{\code{\link[qtl]{argmax.geno}},
  \code{\link[qtl]{find.errors}},
  \code{\link[qtl]{calc.errorlod}},
  \code{\link[qtl]{top.errorlod}},
  \code{\link[qtl]{plot.errorlod}} }

\keyword{hplot}
 
