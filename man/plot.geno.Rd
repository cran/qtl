\name{plot.geno}
\alias{plot.geno}

\title{Plot observed genotypes, flagging likely errors}

\description{
  Plot the genotypes on a particular chromosome for a set of
  individuals, flagging likely errors.
}

\usage{
plot.geno(x, chr, ind, horizontal=FALSE, cutoff=3.5,
          min.sep=2, cex=1.2, \dots) 
}			    

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link[qtl]{read.cross}} for details. }
  \item{chr}{Chromosome number to plot.}
  \item{ind}{Vector of individuals to plot.  If missing, all individuals
    are plotted.}
  \item{horizontal}{If TRUE, chromosomes are plotted horizontally.}
  \item{cutoff}{Genotypes with error LOD scores above this value are
    flagged as possible errors.}
  \item{min.sep}{Markers separated by less than this value (as a
    percent of the chromosome length) are pulled apart, so that they may
    be distinguished in the picture.}
  \item{cex}{Character expansion for the size of points in the plot.
    Larger numbers give larger points; see \code{\link[qtl]{par}}.}
  \item{\dots}{Ignored at this point.}
}

\details{
  Plots the genotypes for a set of individuals.  Likely errors are
  indicated by red squares.  In a backcross, genotypes AA and AB are
  indicated by white and black circles, respectively.  In an intercross,
  genotypes AA, AB and BB are indicated by white, gray, and black
  circles, respectively, and the partially missing genotypes ``not BB''
  (D in mapmaker) and ``not AA'' (C in mapmaker) are indicated by green
  and orange circles, respectively.
}

\value{None.}

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,19))}
# run calc.genoprob, then calc.errorlod; use same error.prob
hyper <- calc.genoprob(hyper,error.prob=0.01)
hyper <- calc.errorlod(hyper,error.prob=0.01)

# plot the error LOD scores; print those above a specified cutoff 
plot.errorlod(hyper)
plot.errorlod(hyper,chr=1)
top.errorlod(hyper,cutoff=4)

# plot genotype data, flagging genotypes with error LOD > cutoff
plot.geno(hyper, chr=1, ind=160:200, cutoff=7, min.sep=2)

# Note that if markers are extremely tightly linked,
#     any recombination event will have a high error LOD
plot.geno(hyper, chr=1, ind=101:118, min.sep=2)
hyper$geno[[1]]$map
}

\seealso{ \code{\link[qtl]{calc.errorlod}},
  \code{\link[qtl]{plot.errorlod}}, \code{\link[qtl]{top.errorlod}},
  \code{\link[base]{image}} }


\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{hplot}
