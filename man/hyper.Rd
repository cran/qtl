\name{hyper}
\alias{hyper}
\non_function{}
\title{Data on hypertension}
\description{ 
  This is data from a real experiment on hypertension in the mouse.  The
  phenotype is the blood pressure.  See the references below. Note that,
  for most markers, genotypes are available on only the individuals with
  extreme phenotypes.  At many markers, only recombinant individuals
  were typed.  We thank Bev Paigen and Gary Churchill for allowing us to
  distribute these data with the R/qtl package. 
} 
\usage{data(hyper)}
\format{
  The data is an object of class \code{cross}.  See
  \code{\link[qtl]{read.cross}} for details. There are 250 backcross
  individuals typed at 174 markers (actually one contains only missing
  values), including 4 on the X chromosome, with one phenotype (blood
  pressure). }  

\references{
  Sugiyama et al. (2001) Concordance of murine quantitative trait loci
  for salt-induced hypertension with rat and human loci. Genomics 71:70-77.
}
  
\source{Bev Paigen and Gary Churchill (The Jackson Laboratory, Bar
  Harbor, Maine)
  \url{http://www.jax.org/research/churchill/datasets/qtl/qtlarchive} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(hyper)
}

\seealso{\code{\link[qtl]{fake.bc}}, \code{\link[qtl]{fake.f2}},
  \code{\link[qtl]{fake.4way}}, \code{\link[qtl]{listeria}} } 

\keyword{datasets}
