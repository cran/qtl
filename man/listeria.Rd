\name{listeria}
\alias{listeria}
\non_function{}
\title{Data on Listeria monocytogenes susceptibility}
\description{ 
  This is data from a real experiment on susceptibility to Listeria
  monocytogenes infection in the mouse.  The phenotype is the survival
  time (in hours) following infection.  Mice with phenotype 264.0 hours
  may be considered to have recovered from the infection. See the
  references below.  We thank Victor Boyartchuk and William Dietrich for
  allowing us to distribute these data with the R/qtl package.
} 
\usage{data(listeria)}
\format{
  The data is an object of class \code{cross}.  See
  \code{\link[qtl]{read.cross}} for details.
  There are 120 F2 individuals typed at 133 markers, including 2 on the
  X chromosome, with one phenotype.
}

\references{
  Boyartchuk VL, Broman KW, Mosher RE, D'Orazio SEF, Starnbach MN,
  Dietrich WF (2001) Multigenic control of Listeria monocytogenes
  susceptibility in mice. Nature Genetics 27:259-260 \cr
  
  Broman KW, Boyartchuk VL, Dietrich WF (2000) Mapping time-to-death
  quantitative trait loci in a mouse cross with high survival
  rates. Technical Report MS00-04, Department of Biostatistics, Johns
  Hopkins University, Baltimore, MD
}
  
\source{Victor Boyartchuk and William Dietrich (Department of Genetics,
  Harvard Medical School and Howard Hughes Medical Institute) 
  \url{http://genetics.med.harvard.edu/~bdlab} }

\author{Karl W Broman,
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\examples{
data(listeria)
}

\seealso{\code{\link[qtl]{fake.bc}}, \code{\link[qtl]{fake.f2}},
  \code{\link[qtl]{fake.4way}}, \code{\link[qtl]{hyper}} } 

\keyword{datasets}
