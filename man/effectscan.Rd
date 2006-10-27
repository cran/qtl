\name{effectscan}
\alias{effectscan}

\title{Plot allelic effect across the whole genome}

\description{
  This function is used to plot allelic effect at all markers for
  selected chromosomes. For backcross, there will be only one line,
  representing the additive effect. For F2, there will be two lines,
  representing the additive and dominance effects. 
}

\usage{
effectscan(cross, pheno.col=1, chr, ylim, gap=25,
           col=c("black","blue","red"), lty=c(1,2,3), lwd=2,
           mtick=c("line", "triangle"), main, add.legend=TRUE,
           \dots)
}

\arguments{
  \item{cross}{An object of class \code{cross}.}
  \item{pheno.col}{Column number in the phenotype matrix which to be
    drawn in the plot.}
  \item{chr}{Chromosome(s) to be drawn in the plot (optional).}
  \item{ylim}{y-axis limit.}
  \item{gap}{Gap separating chromosomes (in cM).}
  \item{col}{Line colors.}
  \item{lty}{Line types.}
  \item{lwd}{Line widths.}
  \item{mtick}{Tick mark type for markers.}
  \item{main}{Figure title.}
  \item{add.legend}{A logical value to indicate add legend or not.}
  \item{\dots}{Passed to the function \code{\link[graphics]{plot}} when it
    is called.}
}
\value{None.}

\details{
  The function calls \code{\link[qtl]{effectplot}} repeatedly for all
  markers on the selected chromosomes.  The results of
  \code{\link[qtl]{sim.geno}} are required for taking account of missing
  genotype information.
}

\examples{
data(fake.f2)
# allelic effect on whole genome
effectscan(fake.f2, ylim=c(-6,3))
# on first 4 chromosomes, use triangle as marker ticks
effectscan(fake.f2, chr=1:4, ylim=c(-6,3), mtick="triangle")
}

\author{Hao Wu}

\seealso{ \code{\link[qtl]{effectplot}}, \code{\link[qtl]{plot.pxg}},
  \code{\link[qtl]{sim.geno}} } 

\keyword{hplot}
