\name{plot.map}
\alias{plot.map}

\title{Plot genetic map}

\description{
  Plot genetic map of marker locations for all chromosomes.
}

\usage{plot.map(x, map2, horizontal=FALSE, \dots)}

\arguments{
  \item{x}{A list whose components are vectors of marker locations.  A
    \code{cross} object may be given instead, in which case the genetic
    map it contains is used.}
  \item{map2}{An optional second genetic map with the same number of
    chromosomes and markers as the first.  As with the first argument, a
    \code{cross} object may be given instead.  If this argument is
    given, a comparison of the two genetic maps is plotted.}
  \item{horizontal}{Specifies whether the chromosomes should be
    plotted horizontally.}
  \item{\dots}{Ignored at this point.}
}

\value{None.}

\details{
  Plots the genetic map for each chromosome, or a comparison of the
  genetic maps if two maps are given.

  For a sex-specific map, female and male maps are plotted against one
  another.  For two sex-specific maps, female maps and male maps are
  plotted against one another.
}

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=18:19)}
plot.map(fake.bc)
plot.map(fake.bc,horizontal=TRUE)
newmap <- est.map(fake.bc)
plot(newmap)
plot.map(fake.bc, newmap)
}

\seealso{ \code{\link[qtl]{est.map}}, \code{\link[qtl]{plot.cross}} }

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{hplot}
