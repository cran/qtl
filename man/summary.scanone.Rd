\name{summary.scanone}
\alias{summary.scanone}

\title{Summarize the results of one or several one-QTL genome scans}

\description{
  Print the rows of the output from \code{\link[qtl]{scanone}} that
  correspond to the maximum LOD for each chromosome, provided that they
  exceed some specified thresholds.
}

\usage{
summary.scanone(\dots, threshold=0)
}

\arguments{
 \item{\dots}{One or several objects output by the function
   \code{\link[qtl]{scanone}}.}
 \item{threshold}{LOD score thresholds. Only peaks with LOD score above
   this value will be returned. This could be a scalar, a vector or a
   matrix. See below section for details.}
}

\details{
  This function is used to report locus deemed interesting for a one-QTL
  genome scan. Users can provide one or more results from
  \code{\link[qtl]{scanone}} and one or more thresholds corresponding to
  different significance levels. For multiple scanone results case, a loci
  will be reported if any of the LOD scores exceeds the minimum of the
  corresponding threshold.
  
  A column is significance code will be printed for each column of LOD
  score. It will use "stars" to visualize the significance level. That is,
  if the LOD score exceeds the lowest threshold, it will have one star
  (*). If it exceeds the second lowest threshold, it will have two stars
  (**) and so on. 
  
  The input \code{threshold} could be a scalar, a vector or a matrix. The
  way \code{threshold} works is like the following: 
  \enumerate{
    \item One \code{\link[qtl]{scanone}} result case.
    \itemize{
      \item Threshold is a scalar. This is the easiest case. Locus with
      LOD scores exceed this value will be reported.
      \item Threshold is a vector. Different values correspond to
      different significance level.
    }
    \item Multiple \code{\link[qtl]{scanone}} result case.
    \itemize{
      \item Threshold is a scalar. The code will assume the same threshold
      value will apply to all scan results and there will be only one
      significance level.
      \item Threshold is a vector. The code will assume that each value in
      the vector corresponds to a scanone result, and there will be only
      one significance level. Note that the length of the vector must be
      the same as the number of scanone results.
      \item Threshold is a matrix. Each row of the matrix corresponds to a
      significance level and each column corresponds to a scan result.
    }
  }
  
  Usually the threshold should be obtained from the permutation
  test. Users can use \code{\link[stats]{quantile}} to get a list of
  threhold values for different significance levels.
}

\value{
  An object of class \code{summary.scanone}, to be printed by
  \code{print.summary.scanone}.  This is a data.frame with one row per
  chromosome, corresponding to the maximum LOD scores.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu}; Hao Wu,
  \email{hao@jax.org} }

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2)
# scanones on observed data
scan1 <- scanone(fake.f2, pheno.col=1)
# permutation tests
\dontrun{out1 <- scanone(fake.f2, pheno.col=1, n.perm=100)
}\dontshow{out1 <- scanone(fake.f2, pheno.col=1, n.perm=3)
}
###### summary one result #######
# one sig level
summary(scan1, threshold=quantile(out1, 0.95))
# multiple sig level
summary(scan1, threshold=quantile(out1, c(0.63,0.95,0.99)))

# scan the other phenotype
scan2 <- scanone(fake.f2, pheno.col=2)
\dontrun{out2 <- scanone(fake.f2, pheno.col=2, n.perm=100)
}\dontshow{out2 <- scanone(fake.f2, pheno.col=2, n.perm=3)
}

##### summary multiple results
# one sig level
summary(scan1,scan2, threshold=2)
# multiple sig level, same threshold
summary(scan1,
        scan2,threshold=matrix(quantile(out1,c(0.37,0.95,0.99)),ncol=1))   
# one sig level, different threshold
summary(scan1, scan2,threshold=c(3,0.8))
# multiple sig level, different threshold
threshold <- matrix(c(quantile(out1,probs=c(0.37,0.95,0.99)),0.5,0.8,1.2), ncol=2) 
summary(scan1, scan2, threshold=threshold)
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}},
          \code{\link[qtl]{max.scanone}}, \code{\link[stats]{quantile}}}

\keyword{print}
