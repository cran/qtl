\name{effectplot}
\alias{effectplot}

\title{Plot phenotype means against genotypes at one or two markers.}

\description{
  Plot the phenotype means for each group defined by the genotypes at 
  one or two markers (or the values at a discrete covariate).
}

\usage{
effectplot(cross, pheno.col=1, mname1, mark1, geno1, mname2, mark2,
           geno2, main, ylim, add.legend=TRUE, draw=TRUE,
           var.flag=c("pooled","group"))
}

\arguments{
  \item{cross}{An object of class \code{cross}.}
  \item{pheno.col}{Column number in the phenotype matrix which to be
    drawn in the plot.}
  \item{mname1}{Name for the first marker.}
  \item{mark1}{Genotype data for the first marker.  If unspecified,
    genotypes will be taken from the data in the input cross object,
    using the name specified in \code{mname1}.}
  \item{geno1}{Names for the genotypes (or classes in a covariate).}
  \item{mname2}{Name for the second marker (optional).}
  \item{mark2}{Like \code{mark1} (optional).}
  \item{geno2}{Names for the genotypes (or classes in a covariate).}
  \item{main}{Optional figure title.}
  \item{ylim}{Y-axis limits.}
  \item{add.legend}{A logical value to indicate add legend or not.}
  \item{draw}{A logical value to indicate generate the plot or not. If
    FALSE, no figure will be plot and this function can be used to
    calculate the group means and standard errors.}
  \item{var.flag}{The method to calculate the group variance. "pooled"
    means to use the pooled variance and "group" means to calculate from
    individual group.}
}

\value{
  A data.frame containing the phenotype means and standard deviations
  for each group.
}

\details{
  In the plot, the Y-axis is the phenotype. X-axis is the genotype for
  the first marker.  In the case of one marker, there will be one line
  in the figure. In the case of two markers, the x-axis is for different
  genotypes of the second marker, and the first marker is represented by
  lines in different colors.  Error bars are plotted at \eqn{\pm}{+/-} 1
  SE. 
}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(1, 13, "X"))}

# impute gneotype data
\dontrun{fake.f2.impute <- sim.geno(fake.f2, step=5, n.draws=64)
}\dontshow{fake.f2.impute <- sim.geno(fake.f2, step=5, n.draws=8)
}

########################################
# one marker plots
########################################
### plot of genotype-specific phenotype means for 1 marker
mname <- find.marker(fake.f2, 1, 37) # marker D1M437
effectplot(fake.f2, pheno.col=1, mname1=mname)

### plot a pseudomarker
effectplot(fake.f2.impute, pheno.col=1, mname1=mname)

### plot a phenotype
# Plot of sex-specific phenotype means,
# note that "sex" must be a phenotype name here
effectplot(fake.f2, mname1="sex", geno1=c("F","M"))
# alternatively I can do it this way
sex <- fake.f2$pheno[,2]
effectplot(fake.f2, mname1="Sex", mark1=sex, geno1=c("F","M"))

########################################
# two markers plots
########################################

### plot two markers
# plot of genotype-specific phenotype means for 2 markers
mname1 <- find.marker(fake.f2, 1, 37) # marker D1M437
mname2 <- find.marker(fake.f2, 13, 24) # marker D13M254
effectplot(fake.f2, mname1=mname1, mname2=mname2)

### plot two pseudomarkers
pmname1 <- "c1.loc35"
pmname2 <- "c13.loc25"
effectplot(fake.f2.impute, mname1=pmname1, mname2=pmname2)

### plot one marker and one imputed pseudomarker
mark1 <- fake.f2.impute$geno[[1]]$draws[,11,]
effectplot(fake.f2, mname1="D1M437", mark1=mark1, mname2="D13M254")

### Plot of sex- and genotype-specific phenotype means 
mname <- find.marker(fake.f2, 13, 24) # marker D13M254
# sex and a marker
effectplot(fake.f2, mname1=mname, mname2="Sex",
           mark2=fake.f2$pheno$sex, geno2=c("F","M"))
# sex and a pseudomarker
effectplot(fake.f2.impute, mname1=mname, mname2="Sex",
           mark2=fake.f2$pheno$sex, geno2=c("F","M"))

# Same as above, switch role of sex and the marker
# sex and marker
effectplot(fake.f2, mname1="Sex", mark1=fake.f2$pheno$sex,
           geno1=c("F","M"), mname2=mname)
# sex and pseudomarker
effectplot(fake.f2.impute, mname1="Sex", mark1=fake.f2$pheno$sex,
           geno1=c("F","M"), mname2=mname)

# X chromosome marker
mname <- find.marker(fake.f2, "X", 14) # marker DXM66
effectplot(fake.f2, mname1=mname)

# Two markers, including one on the X
mname1 <- find.marker(fake.f2, 13, 24) # marker D13M254
mname2 <- find.marker(fake.f2, "X", 14) # marker DXM66
effectplot(fake.f2, mname1=mname1, mname2=mname2)
}

\author{Hao Wu, \email{hao@jax.org}; Karl W Broman, \email{kbroman@jhsph.edu} }

\seealso{ \code{\link[qtl]{plot.pxg}}, \code{\link[qtl]{find.marker}},
  \code{\link[qtl]{effectscan}}} 

\keyword{hplot}
