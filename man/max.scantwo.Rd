\name{max.scantwo}
\alias{max.scantwo}

\title{Maximum peak in two-dimensional genome scan}

\description{
  Print the maximum joint LOD and interaction LODs in the results of
  \code{scantwo}. 
}

\usage{
max.scantwo(\dots, na.rm=TRUE)
}

\arguments{
 \item{\dots}{An object of class \code{scantwo}, the output of the
   function \code{\link[qtl]{scantwo}}.  Note that while multiple
   objects may be input, only the first is considered.}
 \item{na.rm}{A logical indicating whether missing values should be
   removed.}
}

\value{
  An object of class \code{summary.scantwo}, to be printed by
  \code{print.summary.scantwo}.  Pairs of loci that have the maximum
  joint and interaction LOD, in the two-dimensional scan.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} } 

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2, step=5)
out.2dim <- scantwo(fake.f2, method="hk", chr=c(1,6,13))
max(out.2dim)
}

\seealso{ \code{\link[qtl]{scantwo}}, \code{\link[qtl]{plot.scantwo}},
          \code{\link[qtl]{summary.scantwo}} }

\keyword{print}
