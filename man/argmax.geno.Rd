\name{argmax.geno}
\alias{argmax.geno}

\title{Reconstruct underlying genotypes}

\description{
  Uses the Viterbi algorithm to identify the most likely sequence
  of underlying genotypes, given the observed multipoint marker data,
  using a model for genotyping errors. 
}

\usage{
argmax.geno(cross, step=0, off.end=0, error.prob=0, 
            map.function=c("haldane","kosambi","c-f"))
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{step}{Maximum number of cM between positions at which the
 genotypes are to be reconstructed, though for \code{step = 0},
 genotypes are reconstructed only at the marker locations.}
 \item{off.end}{Distance (in cM) at which to carry the genotype
 reconstructions past the p and q terminal markers on each
 chromosome.} 
 \item{error.prob}{Assumed genotyping error rate used in the calculation
 of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
 map function when converting genetic distances into recombination
 fractions.} 
}

\details{
  We use the Viterbi algorithm to calculate
  \eqn{\arg \max_v \Pr(g = v | O)}{arg max_v Pr(g = v | O)} where
  \eqn{g} is the underlying sequence of genotypes and \eqn{O} is the
  observed marker genotypes.

  This done by calculating
  \eqn{\gamma_k(v_k) = \max_{v_1, \ldots, v_{k-1}} \Pr(g_1 = v_1,
    \ldots, g_k = v_k, O_1, \ldots, O_k)}{% 
    gamma_k(v_k) = max{v_1, \ldots, v_{k-1}} Pr(g_1 = v_1,
    \ldots, g_k = v_k, O_1, \ldots, O_k)}
  for \eqn{k = 1, \ldots, n} and then tracing back through the sequence.

  Calculations are done within the C function \code{argmax_geno}.

  Attributes "error.prob", "step", and "off.end" are set to the values of
  the corresponding arguments, for later reference.

  \emph{Note}: the results of \code{argmax.geno()} depend greatly on the
  choice of the value of the argument \code{step}.  This is a sad truth
  of the result of the Viterbi algoriothm.  Further note that if several
  sequences are all most likely, our method of choosing a random such is
  flawed.  
}

\value{
  The \code{cross} object in the input is returned with the reconstructed
  genotypes added.  Recall that the \code{cross$geno} component is a
  list whose elements correspond to chromosomes and which are themselves
  lists with components \code{data} and \code{map}.  For each
  chromosomes, an additional component, \code{argmax}, is added.  This
  is a matrix of size [n.ind x n.pos], where n.pos is the number of
  positions at which the reconstructed genotypes were desired,
  containing the most likely sequences of underlying genotypes.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2 <- argmax.geno(fake.f2,step=2,off.end=5)
\testonly{data(fake.bc)}
\testonly{fake.bc <- argmax.geno(fake.bc,step=0,off.end=0)}
}

\references{
  K Lange (1999) Numerical analysis for statisticians
  Springer-Verlag, New York. Sec 23.3.

  LR Rabiner (1989) A tutorial on hidden Markov models and selected
  applications in speech recognition.  Proceedings of the IEEE
  77:257-286.
}
  

\seealso{\code{\link[qtl]{sim.geno}}, \code{\link[qtl]{calc.genoprob}}}

\keyword{utilities}
