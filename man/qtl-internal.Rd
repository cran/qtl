\name{qtl-internal}
\alias{calc.pairprob}
\alias{checkcovar}
\alias{convert.cross}
\alias{create.map}
\alias{discan}
\alias{imf.cf}
\alias{imf.h}
\alias{imf.k}
\alias{mf.cf}
\alias{mf.h}
\alias{mf.k}
\alias{pull.chr}
\alias{read.cross.csv}
\alias{read.cross.gary}
\alias{read.cross.karl}
\alias{read.cross.mm}
\alias{ripple.perm1}
\alias{ripple.perm2}
\alias{ripple.perm.sub}
\alias{scanone.perm}
\alias{scantwo.perm}
\alias{sim.cross.4way}
\alias{sim.cross.bc}
\alias{sim.cross.f2}
\alias{print.summary.cross}
\alias{print.summary.ripple}
\alias{print.summary.scanone}
\alias{print.summary.scantwo}
\alias{vbscan}
\alias{write.cross.csv}
\alias{write.cross.mm}

\title{Internal qtl functions}

\description{
  Internal qtl functions.
}

\usage{
calc.pairprob(cross, step=0, off.end=0, error.prob=0, 
              map.function=c("haldane","kosambi","c-f"))
checkcovar(cross, pheno.col, addcov, intcov)
convert.cross(cross)
create.map(map, step, off.end)
discan(cross, pheno.col=1, method=c("em","mr"),
       maxit=4000, tol=1e-4)
imf.cf(r)
imf.h(r)
imf.k(r)
mf.cf(d)
mf.h(d)
mf.k(d)
pull.chr(cross, chr)
read.cross.csv(dir, file, sep=",", na.strings="-",
               genotypes=c("A","H","B","C","D"),
               estimate.map=TRUE)
read.cross.gary(dir, genfile, mnamesfile, chridfile, phefile,
                pnamesfile, mapfile)
read.cross.karl(dir, genfile, mapfile, phefile)
read.cross.mm(dir, rawfile, mapfile, estimate.map)
ripple.perm1(n)
ripple.perm2(n)
ripple.perm.sub(x, mat)
scanone.perm(cross, pheno.col=1, model=c("normal","binary","2part","np"),
             method=c("em","imp","hk","mr"), addcov=NULL, intcov=NULL,
             upper=TRUE, ties.random=FALSE, start=NULL, maxit=4000,
             tol=1e-4, n.perm=1000, trace=TRUE)
scantwo.perm(cross, pheno.col=1, method=c("em","imp","hk","mr"),
             addcov=NULL, intcov=NULL, incl.markers=FALSE,
             maxit=4000, tol=1e-4, trace=FALSE, n.perm=1000)
sim.cross.4way(map, model, n.ind, error.prob, missing.prob,
               partial.missing.prob, keep.errorind, map.function)
sim.cross.bc(map, model, n.ind, error.prob, missing.prob,
             keep.errorind, map.function)
sim.cross.f2(map, model, n.ind, error.prob, missing.prob,
             partial.missing.prob, keep.errorind, map.function)
print.summary.cross(x, \dots)
print.summary.ripple(x, \dots)
print.summary.scanone(x, \dots)
print.summary.scantwo(x, \dots)
vbscan(cross, pheno.col=1, upper=FALSE, method="em",
       maxit=4000, tol=1e-4)
write.cross.csv(cross, filestem="data", chr, digits=5)
write.cross.mm(cross, filestem="data", chr, digits=5)
}

\details{
  These are generally not to be called by the user.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{internal}
