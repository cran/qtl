\name{listeria}
\alias{listeria}
\non_function{}

\title{Data on Listeria monocytogenes susceptibility}

\description{ 
  Data from an experiment on susceptibility to \emph{Listeria
  monocytogenes} infection in the mouse.  
} 

\usage{data(listeria)}

\format{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}}
  for details. 
}

\details{
  There are 120 F2 individuals typed at 133 markers, including 2 on the
  X chromosome, with one phenotype.

  The phenotype is the survival time (in hours) following infection.
  Mice with phenotype 264 hours may be considered to have recovered
  from the infection. See the references below.  
}

\references{
  Boyartchuk, V. L., Broman, K. W., Mosher, R. E., D'Orazio
  S. E. F., Starnbach, M. N. and Dietrich, W. F. (2001) Multigenic
  control of \emph{Listeria monocytogenes} susceptibility in
  mice. \emph{Nature Genetics} \bold{27}, 259--260.
  
  Broman, K. W., Boyartchuk, V. L. and Dietrich, W. F. (2000) Mapping
  time-to-death quantitative trait loci in a mouse cross with high
  survival rates. Technical Report MS00-04, Department of Biostatistics,
  Johns Hopkins University, Baltimore, MD.
}
  
\source{
  Victor Boyartchuk and William Dietrich (Department of Genetics,
  Harvard Medical School and Howard Hughes Medical Institute)
  \url{http://genetics.med.harvard.edu/~bdlab}
}

\examples{
data(listeria)
\dontshow{listeria <- subset(listeria,chr=c(1,5,6,13,15))}
# Summaries
summary(listeria)
plot(listeria)

# Take log of phenotype
listeria$pheno[,2] <- log2(listeria$pheno[,1])
names(listeria$pheno)[2] <- "logSurv"
plot(listeria)

# Genome scan with a two-part model, using log survival
listeria <- calc.genoprob(listeria, step=2)
out <- scanone(listeria, model="2part", method="em",
               upper=TRUE, pheno.col=2)

# Summary of the results
summary(out, 4)

# Plot LOD curves for interesting chromosomes
#     (The two-part model gives three LOD scores)
plot(out, out, out, chr=c(1,5,6,13,15), lodcolumn=3:5,
     lty=1, col=c("black","red","blue"))
}

\seealso{ \code{\link[qtl]{fake.bc}}, \code{\link[qtl]{fake.f2}},
  \code{\link[qtl]{fake.4way}}, \code{\link[qtl]{hyper}},
  \code{\link[qtl]{bristle3}}, \code{\link[qtl]{bristleX}} }

\keyword{datasets}
