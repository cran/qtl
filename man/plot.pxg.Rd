\name{plot.pxg}
\alias{plot.pxg}

\title{Plot phenotypes versus marker genotypes.}

\description{
  Plot the phenotype values versus the genotypes at a marker.
}

\usage{
plot.pxg(x, marker, pheno.col=1, jitter=1, infer=TRUE, \dots)
}			    

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link[qtl]{read.cross}} for details. }
  \item{marker}{A marker name (a character string).}
  \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.}
  \item{jitter}{A positive number indicating how much to spread out the
    points horizontally.  (Larger numbers correspond to greater spread.)}
  \item{infer}{If TRUE, missing genotypes are filled in with a single
    random imputation and plotted in red; if FALSE, only individuals
    typed at the specified marker are plotted.}
  \item{\dots}{Ignored at this point.}
}

\details{
  Plots the phenotype data against the genotypes at the specified
  marker.  If infer=TRUE, the genotypes of individuals that were not
  typed is inferred based the genotypes at linked markers via a single
  imputation from \code{\link[qtl]{sim.geno}}; these points are plotted
  in red.  For each genotype, the phenotypic mean is plotted, with error
  bars at \eqn{\pm}{+/-} 1 SE.
}

\value{
  Returns a list with two components.  The first component is a nominal p-value
  corresponding to a test of a difference between the genotype-specific
  phenotypic means.  The second component is a matrix with three
  columns: the genotypes at the marker, the phenotypes, and a 0/1
  indicator with 1's corresponding to inferred genotypes.
}

\examples{
data(listeria)
mymarker <- find.marker(listeria, 5, 28) # marker D5M357
plot.pxg(listeria, mymarker)

data(fake.f2)
mname <- find.marker(fake.f2, 1, 37) # marker D1M437
plot.pxg(fake.f2, mname)

mname <- find.marker(fake.f2, "X", 14) # marker DXM66
plot.pxg(fake.f2, mname)
}

\seealso{ \code{\link[qtl]{find.marker}}, \code{\link[qtl]{effectplot}} }


\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{hplot}
