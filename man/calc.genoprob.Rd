\name{calc.genoprob}
\alias{calc.genoprob}

\title{Calculate conditional genotype probabilities}

\description{
  Uses the hidden Markov model technology to calculate the
  probabilities of the true underlying genotypes given the observed
  multipoint marker data, using a model for genotyping errors.
}

\usage{
calc.genoprob(cross, step=0, off.end=0, error.prob=0, 
              map.function=c("haldane","kosambi","c-f"))
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{step}{Maximum number of cM between positions at which the
 genotype probabilities are calculated, though for \code{step = 0},
 calculate probabilites only at the marker locations.}
 \item{off.end}{Distance (in cM) at which to carry the genotype
 probability calculations past the p and q terminal markers on each
 chromosome.} 
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or Carter-Falconer
 map function when converting genetic distances into recombination
 fractions.} 
}

\details{
  Calculations are done within the C function \code{calc_genoprob}. 

  Let \eqn{O_i} denote the observed marker genotype at position \eqn{i},
  and \eqn{g_i} denote the corresponding true underlying genotype.

  We use the forward-backward equations to calculate
  \eqn{alpha_{iv} = \log Pr(O_1, \ldots, O_i, g_i = v)}{%
    alpha[i][v] = log Pr(O_1, \ldots, O_i, g_i = v)}
  and 
  \eqn{beta_{iv} = \log Pr(O_{i+1}, \ldots, O_n | g_i = v)}{%
    beta[i][v] = log Pr(O_{i+1}, \ldots, O_n | g_i = v)}

  We then obtain
  \eqn{Pr(g_i | O_1, \ldots, O_n) = \exp(alpha_{iv} + beta_{iv}) / s}{%
    Pr(g_i | O_1, \ldots, O_n) = exp(alpha[i][v] + beta[i][v]) / s}
  where
  \eqn{s = \sum_v \exp(alpha_{iv} + beta_{iv})}{%
    s = sum_v exp(alpha[i][v] + beta[i][v])}

  In the case of the 4-way cross, with a sex-specific map, we assume a
  constant ratio of female:male recombination rates within the
  inter-marker intervals. 

  Attributes "error.prob", "step", and "off.end" are set to the values of
  the corresponding arguments, for later reference (especially by the
  function \code{\link[qtl]{cacl.errorlod}}).
}

\value{
  The \code{cross} object in the input is returned with the genotype
  probabilities added.  Recall that the \code{cross$geno} component is a
  list whose elements correspond to chromosomes and which are themselves
  lists with components \code{data} and \code{map}.  For each
  chromosomes, an additional component, \code{prob}, is added.  This
  is an array of size [n.ind x n.pos x n.gen] where n.pos is the number
  of positions at which the probabilities were calculated and n.gen = 3
  for an intercross, = 2 for a backcross, and = 4 for a 4-way cross.
}

\references{
  K Lange (1999) Numerical analysis for statisticians
  Springer-Verlag, New York. Sec 23.3.

  LR Rabiner (1989) A tutorial on hidden Markov models and selected
  applications in speech recognition.  Proceedings of the IEEE
  77:257-286.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2,step=2,off.end=5)
\testonly{data(fake.bc)}
\testonly{fake.bc <- calc.genoprob(fake.bc,step=0,off.end=0)}
}

\seealso{\code{\link[qtl]{sim.geno}}, \code{\link[qtl]{argmax.geno}},
         \code{\link[qtl]{calc.errorlod}}}

\keyword{utilities}
