\name{A starting point}
\alias{A starting point}
\title{Introductory comments}
\description{
  A brief introduction to the R/qtl package.
}

\section{New to R?}{
  \itemize{
    \item Be sure to read the R documention and/or tutorials available
    from the R archive (\url{http://cran.r-project.org}).
    \item In order to use the R/qtl package, you must type (within R)
    \code{library(qtl)}.  You may wish to include this in a
    \code{\link[base:Startup]{.First}} function or
    \code{\link[base:Startup]{.Rprofile}} file. 
    \item Use the \code{\link[base]{help.start}} function to start the
    html version of the R help.
    \item Type \code{library(help=qtl)} to get a list of the functions
    in R/qtl.
    \item Use the \code{\link[base]{example}} function to run examples of the 
    various functions in R/qtl.
  }
}


\section{Walk-through of an analysis}{
You can use the \code{\link[base]{library}} function to list available
packages, list functions and data within a package, and load a package.

\code{library()                  # List available packages} \cr
\code{library(help=qtl)          # List functions and data within "qtl"} \cr
\code{library(qtl)               # Load the qtl package}

Use \code{\link[base:help]{?}} or \code{\link[base]{help}} to get help
about a function or dataset, or use \code{\link[base]{help.start}} to
start the html version of R help.

\code{?calc.genoprob} \cr
\code{help(calc.genoprob)} \cr
\code{help.start()}

Use the \code{\link[base]{data}} function to load one of the example
data sets provided with R/qtl, \code{\link[qtl]{hyper}}.

\code{data(hyper)}

The \code{\link[qtl]{hyper}} data set has class \code{cross}.  The
function \code{\link[qtl]{plot.cross}} calls
\code{\link[qtl]{plot.missing}} (to plot a matrix displaying missing
genotypes) and \code{\link[qtl]{plot.map}} (to plot the genetic maps)
and also displays histograms of all phenotypes.  The
\code{\link[qtl]{plot.missing}} function can plot individuals ordered by
their phenotypes; you can see that for most markers, only individuals
with extreme phenotypes were genotyped.

\code{plot(hyper)} \cr
\code{plot.missing(hyper)} \cr
\code{plot.missing(hyper,reorder=TRUE)} \cr
\code{plot.map(hyper)}

To re-estimate the genetic map for an experimental cross, use the
function \code{\link[qtl]{est.map}}.  The function
\code{\link[qtl]{plot.map}}, in addition to plotting a single map, can
plot the comparison of two genetic maps (as long as they are composed of
the same numbers of chromosomes and markers per chromosome).  The
function \code{\link[qtl]{replace.map}} map be used to replace the
genetic map in a cross with a new one.

\code{newmap <- est.map(hyper,error.prob=0.01,print.rf=TRUE)} \cr
\code{plot.map(hyper,newmap)} \cr
\code{hyper <- replace.map(hyper,newmap)}

Before doing QTL analyses, a number of intermediate calculations may
need to be performed.  The function \code{\link[qtl]{calc.genoprob}}
calculates conditional genotype probabilities given the multipoint
marker data.  \code{\link[qtl]{argmax.geno}} calculates the most likely
sequence of underlying genotypes, given the observed marker data.
\code{\link[qtl]{sim.geno}} simulates sequences of genotypes from their
joint distribution, given the observed marker data.

These three functions return a modified version of the input cross, with
the intermediate calculations included.  

\code{hyper <- calc.genoprob(hyper,step=1,off.end=5,error.prob=0.01)} \cr
\code{hyper <- argmax.geno(hyper,error.prob=0.01)} \cr
\code{hyper <- sim.geno(hyper,step=1,off.end=5,n.draws=16,error.prob=0.01)}

The function \code{\link[qtl]{calc.errorlod}} may be used to assist in
identifying possibly genotyping errors; it calculates the Lincoln and
Lander error lod scores.

\code{hyper <- calc.errorlod(hyper)} \cr
\code{plot.errorlod(hyper)} \cr
\code{top.errorlod(hyper,cutoff=3)}

The function \code{\link[qtl]{est.rf}} estimates the recombination
fraction between each pair of markers, and calculates a LOD score
testing \eqn{r = 0.5}.  This is useful for identifying markers that are
placed on the wrong chromosome.  Note that since, for these data, many
markers were typed only on recombinant individuals, the pairwise
recombination fractions show rather odd patterns.

\code{hyper <- est.rf(hyper)} \cr
\code{plot.rf(hyper)} \cr
\code{plot.rf(hyper,c(1,5,13))}

The function \code{\link[qtl]{scanone}} performs a genome scan with a
single QTL model.  With \code{method="anova"}, analysis of variance at
each marker is performed; any individuals with missing genotypes are
discarded. With \code{method="im"}, Lander and Botstein's interval
mapping is perfomed, while with \code{method="hk"}, the Haley-Knott
approximation to interval mapping is performed. These two methods
require the results from \code{\link[qtl]{calc.genoprob}}.

\code{out.anova <- scanone(hyper,method="anova")} \cr
\code{out.im <- scanone(hyper,method="im")}
\code{out.hk <- scanone(hyper,method="hk")}

The output of \code{\link[qtl]{scanone}} is a matrix with class
\code{scanone}.  The function \code{\link[qtl]{plot.scanone}} may be
used to plot the results, and may plot three sets of results against each
other, as long as they conform appropriately.

\code{plot(out.im)} \cr
\code{plot(out.im,out.hk,c(1,4),lty=1,col=c("black","blue","red"))}

The function \code{\link[qtl]{summary.scanone}} may be used to list
information on the peak LOD for each chromosome

\code{summary(out.im,1.7)}
}


\note{
  Also see the package BQTL (Bayesian QTL mapping toolkit) by
  Charles C. Berry.  His package gave me the idea to have this "A
  starting point" help page.
}

  
\author{Karl W Broman, 
\email{kbroman@jhsph.edu} \cr
\url{http://biosun01.biostat.jhsph.edu/~kbroman/software/qtl.html} }

\keyword{models}
