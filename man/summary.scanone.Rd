\name{summary.scanone}
\alias{summary.scanone}
\title{Print summary of the scanone output}
\description{
  Print the rows of the output from the function \code{scanone} that
  correspond to the maximum LOD for each chromosome. 
}

\usage{
summary.scanone(object,threshold=0,\dots)
}

\arguments{
 \item{object}{An object of class \code{scanone}, the output of the
   function \code{\link[qtl]{scanone}}.  This is a data.frame whose
   first column is the chromosome number and third column is the LOD
   score.}
 \item{threshold}{Only peaks with LOD score above this value will be returned.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.scanone}, to be printed by
  \code{print.summary.scanone}.  This is a data.frame with one row per
  chromosome, that at which the LOD score is a maximum.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2,step=1,off.end=5)
output1 <- scanone(fake.f2,method="im")
summary(output1)
summary(output1,3)
}

\seealso{\code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}},
         \code{\link[qtl]{vbscan}} }

\keyword{print}
