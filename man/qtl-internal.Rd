\name{qtl-internal}
\alias{addqtl}
\alias{adjust.rf.ri}
\alias{calc.pairprob}
\alias{checkcovar}
\alias{convert.cross}
\alias{create.map}
\alias{discan}
\alias{dropqtl}
\alias{fixXdata}
\alias{getsex}
\alias{getgenonames}
\alias{imf.cf}
\alias{imf.h}
\alias{imf.k}
\alias{imf.m}
\alias{locatemarker}
\alias{makeSSmap}
\alias{mf.cf}
\alias{mf.h}
\alias{mf.k}
\alias{mf.m}
\alias{parseformula}
\alias{pull.chr}
\alias{read.cro.qtlcart}
\alias{read.cross.csv}
\alias{read.cross.gary}
\alias{read.cross.karl}
\alias{read.cross.mm}
\alias{read.cross.qtlcart}
\alias{read.cross.qtx}
\alias{read.map.qtlcart}
\alias{read.maps.mm}
\alias{ripple.perm1}
\alias{ripple.perm2}
\alias{ripple.perm.sub}
\alias{scanone.perm}
\alias{scantwo.perm}
\alias{sim.cross.4way}
\alias{sim.cross.bc}
\alias{sim.cross.f2}
\alias{print.summary.cross}
\alias{print.summary.fitqtl}
\alias{print.summary.ripple}
\alias{print.summary.scanone}
\alias{print.summary.scantwo}
\alias{replaceqtl}
\alias{vbscan}
\alias{write.cross.csv}
\alias{write.cross.gary}
\alias{write.cross.mm}
\alias{write.cross.qtlcart}

\title{Internal qtl functions}

\description{
  Internal qtl functions.
}

\usage{
addqtl(cross, qtl, add.chr, add.pos, add.name, map)
adjust.rf.ri(r, type=c("self","sib"), expand=TRUE)
calc.pairprob(cross, step=0, off.end=0, error.prob=0, 
              map.function=c("haldane","kosambi","c-f","morgan"))
checkcovar(cross, pheno.col, addcovar, intcovar)
convert.cross(cross)
create.map(map, step, off.end)
discan(cross, pheno.col=1, method=c("em","mr"),
       x.treatment, maxit=4000, tol=1e-4)
dropqtl(qtl, drop)
fixXdata(type=c("f2ss","f2","bc"), x.treatment=c("simple","full"),
         sexpgm, geno, prob, draws, pairprob)
getsex(cross)
getgenonames(type=c("f2","bc","f2ss","riself","risib","4way"),
             chrtype=c("A","X"), x.treatment=c("simple","full"),
             sexpgm)
imf.cf(r)
imf.h(r)
imf.k(r)
imf.m(r)
locatemarker(map, pos, chr, flag)
makeSSmap(cross)
mf.cf(d)
mf.h(d)
mf.k(d)
mf.m(d)
parseformula(formula, qtl.dimname, covar.dimname)
pull.chr(cross, chr)
read.cro.qtlcart(file)
read.cross.csv(dir, file, sep=",", na.strings=c("-","NA"),
               genotypes=c("A","H","B","D","C"),
               estimate.map=TRUE)
read.cross.gary(dir, genfile, mnamesfile, chridfile, phefile,
                pnamesfile, mapfile,estimate.map)
read.cross.karl(dir, genfile, mapfile, phefile)
read.cross.mm(dir, rawfile, mapfile, estimate.map=TRUE)
read.cross.qtlcart(dir, crofile, mapfile)
read.cross.qtx(dir, file, estimate.map=TRUE)
read.map.qtlcart(file)
read.maps.mm(mapsfile)
ripple.perm1(n)
ripple.perm2(n)
ripple.perm.sub(x, mat)
scanone.perm(cross, pheno.col=1, model=c("normal","binary","2part","np"),
             method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             x.treatment=c("simple","full"), upper=FALSE, ties.random=FALSE,
             start=NULL, maxit=4000, tol=1e-4, n.perm=1000, trace=TRUE)
scantwo.perm(cross, pheno.col=1,
             method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             x.treatment=c("simple","full"), incl.markers=FALSE,
             maxit=4000, tol=1e-4, trace=FALSE, n.perm=1000)
sim.cross.4way(map, model, n.ind, error.prob, missing.prob,
               partial.missing.prob, keep.errorind, map.function)
sim.cross.bc(map, model, n.ind, error.prob, missing.prob,
             keep.errorind, map.function)
sim.cross.f2(map, model, n.ind, error.prob, missing.prob,
             partial.missing.prob, keep.errorind, map.function)
print.summary.fitqtl(x, \dots)
print.summary.cross(x, \dots)
print.summary.ripple(x, \dots)
print.summary.scanone(x, \dots)
print.summary.scantwo(x, \dots)
replaceqtl(cross, qtl, replace, by.chr, by.pos, by.name, map)
vbscan(cross, pheno.col=1, x.treatment, upper=FALSE,
       method="em", maxit=4000, tol=1e-4)
write.cross.csv(cross, filestem="data", chr, digits=5)
write.cross.gary(cross, chr, digits)
write.cross.mm(cross, filestem="data", chr, digits=5)
write.cross.qtlcart(cross, filestem="data", chr)
}

\details{
  These are generally not to be called by the user.
}

\author{Karl W Broman, \email{kbroman@jhsph.edu} }

\keyword{internal}
