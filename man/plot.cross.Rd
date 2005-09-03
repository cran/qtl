\name{plot.cross}
\alias{plot.cross}

\title{Plot various features of a cross object}

\description{
  Plots grid of the missing genotypes, genetic map, and histograms or
  barplots of phenotypes for the data from an experimental cross.
}

\usage{
\method{plot}{cross}(x, auto.layout=TRUE, pheno, \dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{auto.layout}{If TRUE, \code{par(mfrow)} is set so that all plots
   fit within one figure.}
 \item{pheno}{Vector of numbers or character strings corresponding to
   phenotypes that should be plotted.} 
 \item{\dots}{Ignored at this point.}
}

\value{None.}

\details{
  Calls \code{\link[qtl]{plot.missing}}, \code{\link[qtl]{plot.map}} and
  \code{\link[base]{hist}} or \code{\link[base]{barplot}} to plot the
  missing genotypes, genetic map, and histograms or barplots of all
  phenotypes. 

  If \code{auto.format=TRUE}, \code{par(mfrow)} is used with
  \code{ceiling(sqrt(n.phe+2))} rows and the minimum number of columns
  so that all plots fit on the plotting device.

  Numeric phenotypes are displayed as histograms with approximately
  \eqn{\sqrt{n}+5}{sqrt(n)+5} bins.  Phenotypes that are factors are
  displayed as barplots by calling \code{barplot(table(phe))}.
}    

\examples{
data(fake.bc)
plot(fake.bc)
}

\seealso{ \code{\link[qtl]{plot.missing}}, \code{\link[qtl]{plot.map}}, 
  \code{\link[base]{hist}}, \code{\link[base]{barplot}},
  \code{\link[base]{table}} }

\author{Karl W Broman, \email{kbroman@jhsph.edu}; Brian Yandell }

\keyword{hplot}
