\name{fake.bc}
\alias{fake.bc}
\non_function{}

\title{Simulated data for a backcross}

\description{ 
  Simulated data for a backcross, obtained using
  \code{\link[qtl]{sim.cross}}.   
} 

\usage{data(fake.bc)}

\format{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}}
  for details. 
}


\details{
  There are 400 backcross individuals typed at 91 markers and with two
  phenotype and two covariates (sex and age).  

  The two phenotype are due to four QTLs, with no epistasis.  There is
  one on chromosome 2 (at 30 cM), two on chromosome 5 (at 10 and 50 cM),
  and one on chromosome 10 (at 30 cM).  The QTL on chromosome 2 has an
  effect only in the males (sex=1); the two QTLs on chromosome 5 have
  effect in coupling for the first phenotype and in repulsion for the
  second phenotype. Age has an effect of increasing the phenotypes.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc,chr=c(2,5,10))}
summary(fake.bc)
plot(fake.bc)

# genome scans without covariates
fake.bc <- calc.genoprob(fake.bc, step=2.5)
out1.nocovar <- scanone(fake.bc, pheno.col=1)
out2.nocovar <- scanone(fake.bc, pheno.col=2)

# genome scans with covariates
ac <- fake.bc$pheno[,c("sex","age")]
ic <- fake.bc$pheno[,"sex"]
out1.covar <- scanone(fake.bc, pheno.col=1,
                      addcovar=ac, intcovar=ic)
out2.covar <- scanone(fake.bc, pheno.col=2,
                      addcovar=ac, intcovar=ic)

# summaries
summary(out1.nocovar, 3)
summary(out1.covar, 3)
summary(out2.nocovar, 3)
summary(out2.covar, 3)

# plots
plot(out1.covar,out1.nocovar,chr=c(2,5,10))
plot(out2.covar,out2.nocovar,chr=c(2,5,10))
}

\seealso{ \code{\link[qtl]{sim.cross}}, \code{\link[qtl]{fake.4way}},
  \code{\link[qtl]{fake.f2}}, \code{\link[qtl]{listeria}},
  \code{\link[qtl]{hyper}},
  \code{\link[qtl]{bristle3}}, \code{\link[qtl]{bristleX}} }

\keyword{datasets}
