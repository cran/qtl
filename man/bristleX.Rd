\name{bristleX}
\alias{bristleX}
\non_function{}

\title{Data on bristle number in Drosophila}

\description{ 
  Data from bristle number in chromosome X recombinant isogenic lines of
  \emph{Drosophila melanogaster} (fruit flies).
} 

\usage{data(bristleX)}

\format{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}}
  for details. 
}

\details{
  There are 92 chromosome X recombinant isogenic lines, derived from
  inbred lines that were selected for low (A) and high (B) abdominal
  bristle numbers.  A recombinant chromosome X was placed in an isogenic
  low background.

  There are eight phenotypes: the average and SD of the number of
  abdominal and sternopleural bristles in males and females for each
  line.

  Each line is typed at 17 genetic markers on chromosome 3.
}

\references{
  Long, A. D., Mullaney, S. L., Reid, L. A., Fry, J. D., Langley,
  C. H. and MacKay, T. F. C. (1995) High resolution mapping of genetic
  factors affecting abdominal bristle number in \emph{Drosophila
    melanogaster}.  \emph{Genetics} \bold{139}, 1273--1291.
}
  
%\source{
%}

\examples{
data(bristleX)
# Summaries
summary(bristleX)
plot(bristleX)

# genome scan for each of the average phenotypes
bristleX <- calc.genoprob(bristleX, step=2)
out <- vector("list", 4)
for(i in 1:4) out[[i]] <- scanone(bristleX, pheno.col=2*i-1)

# Plot the results
    # maximum LOD score among four phenotypes    
ym <- max(sapply(out, function(a) max(a[,3])))
plot(out[[1]], ylim=c(0,ym))
for(i in 2:4) plot(out[[i]], add=TRUE, col=c("blue","red","green")[i-1])
}

\seealso{ \code{\link[qtl]{bristleX}}, \code{\link[qtl]{listeria}},
  \code{\link[qtl]{fake.bc}}, \code{\link[qtl]{fake.f2}},
  \code{\link[qtl]{fake.4way}}, \code{\link[qtl]{hyper}} }

\keyword{datasets}
