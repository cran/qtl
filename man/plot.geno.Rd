\name{plot.geno}
\alias{plot.geno}
\title{Plot genotypes with indication of likely errors}

\description{
  Plot genotypes on a particular chromosome for a set of individuals,
  with an indication of likely errors.}
}

\usage{
plot.geno(x, chr, ind, horizontal=FALSE, cutoff=2,
         method=c("lod","argmax"), min.sep=1,\dots) 
}			    

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link[qtl]{read.cross}} for details. }
  \item{chr}{Chromosome number to plot.}
  \item{ind}{Vector of individuals to plot.  If missing, all individuals
    are plotted.}
  \item{horizontal}{If TRUE, chromosomes are plotted horizontally.}
  \item{cutoff}{Cutoff on error LOD score (if that method is used);
    genotypes with error LOD scores above this cutoff are flagged as
    errors.}
  \item{method}{Indicates whether to use the error LOD scores calculated
    by \code{\link[qtl]{calc.errorlod}} or the results of
    \code{\link[qtl]{argmax.geno}} and \code{\link[qtl]{find.errors}}.}
  \item{min.sep}{Markers separated by less than this value (as a
    percent of the chromosome length) are pulled apart, so that they may
    be distinguished in the picture.} 
  \item{\dots}{Ignored at this point.}
}

\value{
  A plot of the genotypes for a set of individuals.  Likely errors are
  indicated by red squares.  In a backcross,
  genotypes AA and AB are indicated by white and black circles,
  respectively.  In an intercross, genotypes AA, AB and BB are indicated
  by white, gray, and black circles, respectively, and the partially
  missing genotypes "not BB" (D in mapmaker) and "not AA" (C in
  mapmaker) are indicated by green and orange circles, respectively. 
}

\examples{
data(listeria)
listeria <- calc.genoprob(listeria,error.prob=0.01)
listeria <- calc.errorlod(listeria,error.prob=0.01)
plot.geno(listeria,13)
}

\seealso{\code{\link[base]{image}}, \code{\link[qtl]{plot.cross}} }


\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
