\name{plot.scanone}
\alias{plot.scanone}
\title{Plot output for a single QTL scan}
\description{
  Plot the output from a single QTL scan.  One may specify which
  chromosomes to plot
}
\usage{
plot.scanone(x,output2,output3,chr,incl.markers=TRUE,ylim,
             lty=c(1,2,3),col="black",lwd=2,add=FALSE,gap=25,\dots)
}			    

\arguments{
  \item{x}{A data.frame with at least three columns.  The first
    three columns should be the following:}
  \tabular{ll}{
     \code{chr} \tab Chromosome number \cr
     \code{pos} \tab Position on genetic map \cr
     \code{lod}  \tab LOD score (or other statistic)}
  \item{output2}{Data for another genome scan, like \code{output} (optional).} 
  \item{output3}{Data for another genome scan, like \code{output} (optional).} 
  \item{chr}{Vector specifying which chromosomes to plot.}
  \item{incl.markers}{Indicate whether to plot line segments at the
  marker locations.}
  \item{ylim}{Limits for y-axis [optional].}
  \item{lty}{Line type for output1, 2 and 3; this should be a vector of
    length 1 or 3.}
  \item{col}{Line color for output1, 2 and 3; this should be a vector of
    length 1 or 3.}
  \item{lwd}{Line width for output1, 2 and 3; this should be a vector of
    length 1 or 3.}
  \item{add}{If TRUE, add to a current plot.}
  \item{gap}{Gap separating chromosomes.}
  \item{\dots}{Ignored at this point.}
}
\value{
  Plot of the lod score against marker position for the genome scan.
}

\examples{
data(fake.f2)
fake.f2 <- argmax.geno(fake.f2)
output <- scanone(fake.f2,method="anova")
output2 <- scanone(pull.chr(fake.f2,1),method="im")
plot(output)
plot(output,output2,chr=1)
}

\seealso{\code{\link[qtl]{scanone}}, \code{\link[qtl]{vbscan}},
         \code{\link[qtl]{summary.scanone}} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{hplot}
