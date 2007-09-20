\name{qtl-internal}
\alias{addqtl}
\alias{adjust.rf.ri}
\alias{calc.genoprob.special}
\alias{calc.pairprob}
\alias{checkcovar}
\alias{clean}
\alias{condense}
\alias{convert}
\alias{create.map}
\alias{discan}
\alias{dropqtl}
\alias{dropXcol}
\alias{effectplot.getmark}
\alias{effectplot.calmeanse}
\alias{fitstahl.estp}
\alias{fitstahl.este}
\alias{fitstahl.estpe}
\alias{fitstahl.estp.sub}
\alias{fitstahl.este.sub}
\alias{fitstahl.estpe.sub}
\alias{fixX4write}
\alias{fixXgeno.bc}
\alias{fixXgeno.f2}
\alias{forwsel}
\alias{genotab.em}
\alias{getid}
\alias{getsex}
\alias{getgenonames}
\alias{grab.arg.names}
\alias{imf.cf}
\alias{imf.h}
\alias{imf.k}
\alias{imf.m}
\alias{imf.stahl}
\alias{locatemarker}
\alias{locate.xo}
\alias{makeSSmap}
\alias{max.scantwocondensed}
\alias{mf.cf}
\alias{mf.h}
\alias{mf.k}
\alias{mf.m}
\alias{mf.stahl}
\alias{parseformula}
\alias{plot.effectscan}
\alias{print.cross}
\alias{print.map}
\alias{print.qtl}
\alias{print.scanoneboot}
\alias{print.scantwo}
\alias{print.summary.cross}
\alias{print.summary.fitqtl}
\alias{print.summary.map}
\alias{print.summary.qtl}
\alias{print.summary.ripple}
\alias{print.summary.scanone}
\alias{print.summary.scanoneperm}
\alias{print.summary.scantwo}
\alias{print.summary.scantwo.old}
\alias{print.summary.scantwoperm}
\alias{read.cro.qtlcart}
\alias{read.cross.csv}
\alias{read.cross.csvs}
\alias{read.cross.gary}
\alias{read.cross.karl}
\alias{read.cross.mm}
\alias{read.cross.qtlcart}
\alias{read.cross.qtx}
\alias{read.map.qtlcart}
\alias{read.maps.mm}
\alias{replaceqtl}
\alias{revisecovar}
\alias{revisescantwodf}
\alias{reviseXdata}
\alias{ripple.perm1}
\alias{ripple.perm2}
\alias{ripple.perm.sub}
\alias{scanone.perm}
\alias{scanone.perm.engine}
\alias{scanoneXnull}
\alias{scantwo.perm}
\alias{scantwo.perm.engine}
\alias{sim.bcg}
\alias{sim.cc}
\alias{sim.cross.4way}
\alias{sim.cross.bc}
\alias{sim.cross.f2}
\alias{sim.ril}
\alias{subrousummaryscantwo}
\alias{summary.scantwocondensed}
\alias{vbscan}
\alias{write.cross.csv}
\alias{write.cross.gary}
\alias{write.cross.mm}
\alias{write.cross.qtlcart}

\title{Internal qtl functions}

\description{
  Internal qtl functions. These are generally not to be called by the
  user.
}

\usage{
addqtl(cross, qtl, add.chr, add.pos, add.name, map)
adjust.rf.ri(r, type=c("self","sib"), chrtype=c("A","X"), expand=TRUE)
calc.genoprob.special(cross, error.prob=0.0001,
                      map.function=c("haldane","kosambi","c-f","morgan"))
calc.pairprob(cross, step=0, off.end=0, error.prob=0.0001, 
              map.function=c("haldane","kosambi","c-f","morgan"), map)
checkcovar(cross, pheno.col, addcovar, intcovar, perm.strata, verbose=TRUE)
clean(object)
condense(object)
convert(object)
create.map(map, step, off.end, stepwidth=c("fixed", "variable"))
discan(cross, pheno, method=c("em","mr"),
       addcovar=NULL, intcovar=NULL, maxit=4000, tol=1e-4,
       verbose=FALSE, give.warnings=TRUE)
dropqtl(qtl, drop)
dropXcol(type=c("f2","bc"), sexpgm, cross.attr)
fitstahl.estp(cross, error.prob=0.0001, m=0, tol=1e-4, maxit=4000)
fitstahl.estp.sub(p, cross, error.prob=0.0001, m=0, thetol=1e-4, maxit=4000)
fitstahl.este(cross, m=0, p=0, tol=1e-4, maxit=4000)
fitstahl.este.sub(error.prob, cross, m=0, p=0, thetol=1e-4, maxit=4000)
fitstahl.estpe(cross, m=0, tol=1e-4, maxit=4000)
fitstahl.estpe.sub(x, cross, m=0, thetol=1e-4, maxit=4000)
fixX4write(geno,sex,pgm,crosstype)
fixXgeno.bc(cross)
fixXgeno.f2(cross)
forwsel(x, y, maxsize=7)
genotab.em(dat, tol=1e-6, maxit=10000, verbose=FALSE)
getid(cross)
getsex(cross)
getgenonames(type=c("f2","bc","riself","risib","4way"),
             chrtype=c("A","X"), expandX=c("simple","standard","full"),
             sexpgm, cross.attr)
grab.arg.names(\dots)
imf.cf(r)
imf.h(r)
imf.k(r)
imf.m(r)
imf.stahl(r, m=0, p=0, tol=1e-12, maxit=1000)
locatemarker(map, pos, chr, flag)
locate.xo(cross)
makeSSmap(cross)
\method{max}{scantwocondensed}(object, lodcolumn=1,
    what=c("best", "full", "add", "int"),
    df=FALSE, na.rm=TRUE, \dots)
mf.cf(d)
mf.h(d)
mf.k(d)
mf.m(d)
mf.stahl(d, m=0, p=0)
parseformula(formula, qtl.dimname, covar.dimname)
\method{plot}{effectscan}(x, gap=25, ylim, mtick=c("line","triangle"),
     add.legend=TRUE, alternate.chrid=FALSE, \dots)
\method{print}{cross}(x, \dots)
\method{print}{map}(x, \dots)
\method{print}{qtl}(x, \dots)
\method{print}{scanoneboot}(x, \dots)
\method{print}{scantwo}(x, \dots)
\method{print}{summary.cross}(x, \dots)
\method{print}{summary.fitqtl}(x, \dots)
\method{print}{summary.map}(x, \dots)
\method{print}{summary.qtl}(x, \dots)
\method{print}{summary.ripple}(x, \dots)
\method{print}{summary.scanone}(x, \dots)
\method{print}{summary.scanoneperm}(x, \dots)
\method{print}{summary.scantwo}(x, \dots)
\method{print}{summary.scantwo.old}(x, \dots)
\method{print}{summary.scantwoperm}(x, \dots)
read.cro.qtlcart(file)
read.cross.csv(dir, file, na.strings=c("-","NA"),
               genotypes=c("A","H","B","D","C"),
               estimate.map=TRUE, rotate=FALSE, \dots)
read.cross.csvs(dir, genfile, phefile, na.strings=c("-","NA"),
                genotypes=c("A","H","B","D","C"),
                estimate.map=TRUE, rotate=FALSE, \dots)
read.cross.gary(dir, genfile, mnamesfile, chridfile, phefile,
                pnamesfile, mapfile,estimate.map,na.strings)
read.cross.karl(dir, genfile, mapfile, phefile)
read.cross.mm(dir, rawfile, mapfile, estimate.map=TRUE)
read.cross.qtlcart(dir, crofile, mapfile)
read.cross.qtx(dir, file, estimate.map=TRUE)
read.map.qtlcart(file)
read.maps.mm(mapsfile)
replaceqtl(cross, qtl, replace, by.chr, by.pos, by.name, map)
revisecovar(sexpgm, covar)
revisescantwodf(df)
reviseXdata(type=c("f2","bc"), expandX=c("simple","standard","full"),
            sexpgm, geno, prob, draws, pairprob, cross.attr)
ripple.perm1(n)
ripple.perm2(n)
ripple.perm.sub(x, mat)
scanone.perm(cross, pheno.col=1, model=c("normal","binary","2part","np"),
             method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"), upper=FALSE,
             ties.random=FALSE, start=NULL, maxit=4000,
             tol=1e-4, n.perm=1000, perm.Xsp=FALSE, perm.strata=NULL,
             verbose=TRUE)
scanone.perm.engine(n.perm, cross, pheno.col, model,
                    method, addcovar, intcovar, weights, use,
                    upper, ties.random, start, maxit, tol,
                    verbose, perm.strata)
scanoneXnull(type, sexpgm)
scantwo.perm(cross, pheno.col=1, model=c("normal","binary"),
             method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"),
             incl.markers=FALSE, clean.output=FALSE,
             maxit=4000, tol=1e-4, verbose=FALSE,
             n.perm=1000, perm.strata)
scantwo.perm.engine(n.perm, cross, pheno.col, model,
                    method, addcovar, intcovar, weights, use,
                    incl.markers, clean.output, maxit, tol, verbose,
                    perm.strata)
sim.bcg(n.ind, map, m, p,
        map.function=c("haldane","kosambi","c-f","morgan"))
sim.cc(parents, n.ril=1, type=c("sibmating", "selfing"),
       error.prob=0, missing.prob=0, m=0, p=0, step=0)
sim.cross.4way(map, model, n.ind, error.prob, missing.prob,
               partial.missing.prob, keep.errorind, m, p, map.function)
sim.cross.bc(map, model, n.ind, error.prob, missing.prob,
             keep.errorind, m, p, map.function)
sim.cross.f2(map, model, n.ind, error.prob, missing.prob,
             partial.missing.prob, keep.errorind, m, p, map.function)
sim.ril(map, n.ril=1, type=c("sibmating", "selfing"), n.str=c("2","4","8"),
        m=0, p=0, random.cross=TRUE)
subrousummaryscantwo(object, for.perm=FALSE)
\method{summary}{scantwocondensed}(object, thresholds, 
         what=c("best", "full", "add", "int"),
         perms, alphas, lodcolumn=1, pvalues=FALSE,
         df=FALSE, allpairs=TRUE, \dots)
vbscan(cross, pheno.col=1, upper=FALSE, method="em", maxit=4000,
       tol=1e-4)
write.cross.csv(cross, filestem="data", digits=5, rotate=FALSE,
                split=FALSE)
write.cross.gary(cross, digits)
write.cross.mm(cross, filestem="data", digits=5)
write.cross.qtlcart(cross, filestem="data")
}


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{internal}
