\name{sim.cross}
\alias{sim.cross}
\title{Simulate a QTL experiment}
\description{
  Simulates data for a QTL experiment using a model in which QTLs act additively.
  A backcross, intercross, or phase-known 4-way cross may be simulated.
}
\usage{
sim.cross(map, model=NULL, n.ind=100, type=c("f2", "bc", "4way"),
          error.prob=0, missing.prob=0, partial.missing.prob=0,
          keep.qtlgeno=FALSE,error.ind=TRUE,
	  map.function=c("haldane","kosambi","c-f"))
}
\arguments{
 \item{map}{A list of length n.chr; its components are vectors
 containing the marker locations on each of the chromosomes.}
 \item{model}{A matrix where each row corresponds to a
 different QTL, and gives the chromosome number, cM position and
 effects of the QTL.}
 \item{n.ind}{Number of individuals to simulate.}
 \item{type}{Indicates whether to simulate an intercross ("f2"), a
 backcross ("bc"), or a 4-way cross ("4way").}
 \item{error.prob}{The genotyping error rate to simulate.}
 \item{missing.prob}{The rate of missing genotypes to simulate.}
 \item{partial.missing.prob}{When simulating an intercross or 4-way
   cross, this gives the rate at which markers will be incompletely
   informative (i.e., dominant or recessive).}
 \item{keep.qtlgeno}{If TRUE, genotypes for the simulated QTLs are not
   removed from the marker genotype data or genetic maps.}
 \item{error.ind}{If TRUE, and if \code{error.prob > 0}, the components
   of \code{geno} will each contain a matrix, \code{errors}, with 1's
   indicating a simulated genotyping error.}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
 map function when converting genetic distances into recombination
 fractions.} 
}
\details{
  Meiosis is assumed to exhibit no crossover interference.  If one of
  the chromosomes has class \code{X}, it is assumed to be the X
  chromosome, and is assumed to be segregating in the cross.  Thus, in an
  intercross, it is segregating like a backcross chromosome.  In a 4-way
  cross, a second phenotype, \code{sex}, will be generated.
  In order to ensure that the code is relatively easy to read,
  simulation is actually performed in a set of subfunctions,
  \code{sim.cross.bc}, \code{sim.cross.f2}, and \code{sim.cross.4way}.
}  
\value{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}} for
  details.} 
}
\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\seealso{\code{\link[qtl]{read.cross}}, \code{\link[qtl]{fake.f2}},
         \code{\link[qtl]{fake.bc}} \code{\link[qtl]{fake.4way}} }

\examples{
map <- sim.map()
fake <- sim.cross(map,type="f2", model = rbind(c(1,20,1,2),c(5,20,0.5,0.5)))
}
                    
\keyword{datagen}
