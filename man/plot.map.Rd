\name{plot.map}
\alias{plot.map}
\title{Plot genetic map}
\description{
  Plot genetic map of marker locations for all chromosomes.
}

\usage{
plot.map(x,map2,horizontal=FALSE,\dots)
}			    

\arguments{
  \item{x}{A list of length n.chr, each component being a vector of
    marker locations.  A \code{cross} object may be given instead, in
    which case the genetic map it contains is used.} 
  \item{map2}{An optional second genetic map with the same number of
    chromosomes and markers as the first.  As with \code{map}, a
    \code{cross} object may be given instead.  If this argument is
    given, a comparison of the two genetic maps is plotted.
    \bold{Note}: For sex-specific maps, if two maps are to be plotted,
    we create two plots: the female maps against each other and the male
    maps against each other.} 
  \item{horizontal}{Specifies whether the chromosomes should be
    plotted horizontally.}
  \item{\dots}{Ignored at this point.}
}
\value{
  Plots the genetic map for each chromosome, or a comparison of the
  genetic maps if two maps are given.
}

\examples{
data(fake.bc)
plot.map(fake.bc)
plot.map(pull.map(fake.bc), horizontal=TRUE)
newmap <- est.map(fake.bc)
plot.map(fake.bc,newmap)
plot.map(pull.map(fake.bc),newmap,horizontal=TRUE)
}

\seealso{\code{\link[qtl]{plot.cross}}, \code{\link[qtl]{est.map}}}

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
