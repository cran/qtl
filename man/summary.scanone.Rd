\name{summary.scanone}
\alias{summary.scanone}

\title{Summarize the results of a genome scan}

\description{
  Print the rows of the output from \code{\link[qtl]{scanone}} that
  correspond to the maximum LOD for each chromosome, provided that they
  exceed a specified threshold.
}

\usage{
summary.scanone(object, threshold=0, \dots)
}

\arguments{
 \item{object}{An object of the form output by the function
   \code{\link[qtl]{scanone}}: a data.frame whose first column is the
   chromosome ID and third column is the LOD score.  All other
   columns are ignored.}
 \item{threshold}{Only peaks with LOD score above this value will be
   returned.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.scanone}, to be printed by
  \code{print.summary.scanone}.  This is a data.frame with one row per
  chromosome, corresponding to the maximum LOD scores.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)
\dontshow{listeria <- subset(listeria, chr=18:19)}
listeria <- calc.genoprob(listeria, step=2.5)
out <- scanone(listeria, model="2part")
# All peaks
summary(out)

# Peaks with LOD(p,mu) > 5
summary(out, 5)
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}},
          \code{\link[qtl]{max.scanone}} }

\keyword{print}
