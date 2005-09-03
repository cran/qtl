\name{max.scanone}
\alias{max.scanone}

\title{Maximum peak in genome scan}

\description{
  Print the row of the output from \code{\link[qtl]{scanone}} that
  corresponds to the maximum LOD, genome-wide.
}

\usage{
\method{max}{scanone}(\dots, chr, na.rm=TRUE)
}

\arguments{
 \item{\dots}{An object of the form output by the function
   \code{\link[qtl]{scanone}}: a data.frame whose third column is the
   LOD score.  Note that while multiple objects may be input, only the
   first is considered.}
 \item{chr}{Numeric or character vector indicating chromosomes
   on which to find the maximum LOD score.}
 \item{na.rm}{A logical indicating whether missing values should be
   removed.}
}

\value{
  An object of class \code{summary.scanone}, to be printed by
  \code{print.summary.scanone}.  This is a data.frame with one row,
  corresponding to the maximum LOD peak either genome-wide or for the
  particular chromosome specified.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(listeria)
\dontshow{listeria <- subset(listeria,chr=c(5,13))}
listeria <- calc.genoprob(listeria, step=2.5)
out <- scanone(listeria, model="2part")
# Maximum peak for LOD(p,mu)
max(out)
max(out,chr=13)

# Maximum peak for LOD(p)
max(out[,-3])

# Maximum peak for LOD(mu)
max(out[,-(3:4)])
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}},
          \code{\link[qtl]{summary.scanone}} }

\keyword{print}
