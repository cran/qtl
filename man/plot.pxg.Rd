\name{plot.pxg}
\alias{plot.pxg}

\title{Plot phenotypes versus marker genotypes.}

\description{
  Plot the phenotype values versus the genotypes at a marker or markers.
}

\usage{
plot.pxg(x, marker, pheno.col=1, jitter=1, infer=TRUE,
         pch , ylab, \dots)
}			    

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link[qtl]{read.cross}} for details. }
  \item{marker}{Marker name (a character string; can be a vector).}
  \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.}
  \item{jitter}{A positive number indicating how much to spread out the
    points horizontally.  (Larger numbers correspond to greater spread.)}
  \item{infer}{If TRUE, missing genotypes are filled in with a single
    random imputation and plotted in red; if FALSE, only individuals
    typed at the specified marker are plotted.}
  \item{pch}{Plot symbol.}
  \item{ylab}{Label for y-axis.}
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
  Returns a list with six components.
  \item{Rsq}{R-square, or explained variation.}
  \item{p.value}{Vector of p-values corresponding to an overall test of
    a difference between the genotype-specific phenotypic means, and
    marker main effects and interactions, if there is more than one
    marker. If there are multiple markers and some missing genotypes,
    interactions are collapsed into one component labeled \code{inter}.}
  \item{me}{Means by marker.}
  \item{se}{Standard errors by marker.}
  \item{fit}{The \code{aov} fit.}
  \item{data}{Matrix with at least three columns: the genotypes at each
    marker, the phenotypes, and a 0/1 indicator with 1's corresponding
    to inferred genotypes.}
}

\examples{
data(listeria)
mname <- find.marker(listeria, 5, 28) # marker D5M357
plot.pxg(listeria, mname)

mname2 <- find.marker(listeria, 13, 26) # marker D13Mit147
plot.pxg(listeria, c(mname, mname2))

data(fake.f2)
mname <- find.marker(fake.f2, 1, 37) # marker D1M437
plot.pxg(fake.f2, mname)

mname2 <- find.marker(fake.f2, "X", 14) # marker DXM66
plot.pxg(fake.f2, mname2)

plot.pxg(fake.f2, c(mname,mname2))
}

\seealso{ \code{\link[qtl]{find.marker}}, \code{\link[qtl]{effectplot}},
  \code{\link[qtl]{find.flanking}} }


\author{Karl W Broman, \email{kbroman@jhsph.edu}; Brian Yandell }

\keyword{hplot}
