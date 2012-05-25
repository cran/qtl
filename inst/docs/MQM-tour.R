### R code from vignette source 'MQM-tour.Rnw'

###################################################
### code chunk number 1: setseed
###################################################
set.seed(19696527)


###################################################
### code chunk number 2: MQM-tour.Rnw:144-147
###################################################
library(qtl)
data(map10)
simcross <- sim.cross(map10, type="f2", n.ind=100, missing.prob=0.02)


###################################################
### code chunk number 3: missingdata (eval = FALSE)
###################################################
## geno.image(simcross)


###################################################
### code chunk number 4: missingdataplot
###################################################
geno.image(simcross)


###################################################
### code chunk number 5: MQM-tour.Rnw:221-223
###################################################
# displays warning because MQM ignores the X chromosome in an F2
augmentedcross <- mqmaugment(simcross, minprob=1.0)


###################################################
### code chunk number 6: augment1 (eval = FALSE)
###################################################
## geno.image(augmentedcross)


###################################################
### code chunk number 7: augment1plot
###################################################
geno.image(augmentedcross)


###################################################
### code chunk number 8: MQM-tour.Rnw:249-250
###################################################
augmentedcross <- mqmaugment(simcross, minprob=0.1)


###################################################
### code chunk number 9: augment2 (eval = FALSE)
###################################################
## geno.image(augmentedcross)


###################################################
### code chunk number 10: augment2plot
###################################################
geno.image(augmentedcross)


###################################################
### code chunk number 11: augment3
###################################################
data(multitrait)
msim5 <- simulatemissingdata(multitrait, 5)
msim10 <- simulatemissingdata(multitrait, 10)
msim80 <- simulatemissingdata(multitrait, 80)


###################################################
### code chunk number 12: augment4
###################################################
maug5 <- mqmaugment(msim5)
maug10 <- mqmaugment(msim10, minprob=0.25)
maug80 <- mqmaugment(msim80, minprob=0.80)


###################################################
### code chunk number 13: augmentMinProb
###################################################
maug10minprob <- mqmaugment(msim10, minprob=0.001, verbose=TRUE)
maug10minprobImpute <- mqmaugment(msim10, minprob=0.001, strategy="impute",
                                  verbose=TRUE)
# check how many individuals are expanded:
nind(maug10minprob)
nind(maug10minprobImpute)


###################################################
### code chunk number 14: augment5
###################################################
mqm5 <- mqmscan(maug5)
mqm10 <- mqmscan(maug10)
mqm80 <- mqmscan(maug80)


###################################################
### code chunk number 15: augment5b
###################################################
msim5 <- calc.genoprob(msim5)
one5 <- scanone(msim5)
msim10 <- calc.genoprob(msim10)
one10 <- scanone(msim10)
msim80 <- calc.genoprob(msim80)
one80 <- scanone(msim80)


###################################################
### code chunk number 16: augment6 (eval = FALSE)
###################################################
## op <- par(mfrow = c(2,2))
## plot(mqm5, mqm10, mqm80, col=c("green","blue","red"), main="MQM missing data")
## legend("topleft", c("MQM 5%","MQM 10%","MQM 80%"), col=c("green","blue","red"), lwd=1)
## plot(one5, mqm5, main="5% missing", col=c("black","green"))
## legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)
## plot(one10, mqm10, main="10% missing", col=c("black","blue"))
## legend("topleft", c("scanone","MQM"), col=c("black","blue"), lwd=1)
## plot(one80, mqm80, main="80% missing", col=c("black","red"))
## legend("topleft", c("scanone","MQM"), col=c("black","red"), lwd=1)


###################################################
### code chunk number 17: MQM-tour.Rnw:353-354
###################################################
op <- par(mfrow = c(2,2))
plot(mqm5, mqm10, mqm80, col=c("green","blue","red"), main="MQM missing data")
legend("topleft", c("MQM 5%","MQM 10%","MQM 80%"), col=c("green","blue","red"), lwd=1)
plot(one5, mqm5, main="5% missing", col=c("black","green"))
legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)
plot(one10, mqm10, main="10% missing", col=c("black","blue"))
legend("topleft", c("scanone","MQM"), col=c("black","blue"), lwd=1)
plot(one80, mqm80, main="80% missing", col=c("black","red"))
legend("topleft", c("scanone","MQM"), col=c("black","red"), lwd=1)


###################################################
### code chunk number 18: MQM-tour.Rnw:379-382
###################################################
data(multitrait)
maug_min1 <- mqmaugment(multitrait, minprob=1.0)
mqm_min1 <- mqmscan(maug_min1)


###################################################
### code chunk number 19: MQM-tour.Rnw:389-391
###################################################
mgenop <- calc.genoprob(multitrait, step=5)
m_one <- scanone(mgenop)


###################################################
### code chunk number 20: MQM-tour.Rnw:399-401
###################################################
maug <- mqmaugment(multitrait)
mqm <- mqmscan(maug)


###################################################
### code chunk number 21: MinprobMulti
###################################################
plot(m_one, mqm_min1, col=c("black","green"), lty=1:2)
legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)


###################################################
### code chunk number 22: MQM-tour.Rnw:436-437
###################################################
real_markers <- mqmextractmarkers(mqm)


###################################################
### code chunk number 23: MQM-tour.Rnw:458-463
###################################################
max(mqm)
find.marker(maug, chr=5, pos=35)
multitoset <- find.markerindex(maug, "GH.117C")
setcofactors <- mqmsetcofactors(maug, cofactors=multitoset)
mqm_co1 <- mqmscan(maug, setcofactors)


###################################################
### code chunk number 24: Cofactor4multi (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_co1))
## plot(mqm_co1)


###################################################
### code chunk number 25: MQM-tour.Rnw:480-482
###################################################
# plot after adding first cofactor
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_co1))
plot(mqm_co1)


###################################################
### code chunk number 26: Cofactor4bMULTI (eval = FALSE)
###################################################
## plot(m_one, mqm_co1, col=c("black","green"), lty=1:2)
## legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)


###################################################
### code chunk number 27: MQM-tour.Rnw:499-500
###################################################
plot(m_one, mqm_co1, col=c("black","green"), lty=1:2)
legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)


###################################################
### code chunk number 28: MQM-tour.Rnw:520-524
###################################################
# summary(mqm_co1)
multitoset <- c(multitoset, find.markerindex(maug, find.marker(maug,4,10)))
setcofactors <- mqmsetcofactors(maug,cofactors=multitoset)
mqm_co2 <- mqmscan(maug, setcofactors)


###################################################
### code chunk number 29: twowaycomparison (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_co2))
## plot(mqm_co1, mqm_co2, col=c("blue","green"), lty=1:2)
## legend("topleft", c("one cofactor","two cofactors"), col=c("blue","green"),
##        lwd=1)


###################################################
### code chunk number 30: MQM-tour.Rnw:539-540
###################################################
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_co2))
plot(mqm_co1, mqm_co2, col=c("blue","green"), lty=1:2)
legend("topleft", c("one cofactor","two cofactors"), col=c("blue","green"),
       lwd=1)


###################################################
### code chunk number 31: threewaycomparisonmulti (eval = FALSE)
###################################################
## plot(mqm, mqm_co1, mqm_co2, col=c("green","red","blue"), lty=1:3)
## legend("topleft", c("no cofactors","one cofactor","two cofactors"),
##        col=c("green","red","blue"), lwd=1)


###################################################
### code chunk number 32: MQM-tour.Rnw:559-561
###################################################
# plot closeup of threeway comparison
plot(mqm, mqm_co1, mqm_co2, col=c("green","red","blue"), lty=1:3)
legend("topleft", c("no cofactors","one cofactor","two cofactors"),
       col=c("green","red","blue"), lwd=1)


###################################################
### code chunk number 33: MQM-tour.Rnw:619-623 (eval = FALSE)
###################################################
## autocofactors <- mqmautocofactors(maug, 50)
## mqm_auto <- mqmscan(maug, autocofactors)
## setcofactors <- mqmsetcofactors(maug, 5)
## mqm_backw <- mqmscan(maug, setcofactors)


###################################################
### code chunk number 34: MQM-tour.Rnw:626-630
###################################################
autocofactors <- mqmautocofactors(maug, 50)
mqm_auto <- mqmscan(maug, autocofactors)
setcofactors <- mqmsetcofactors(maug, 5)
mqm_backw <- mqmscan(maug, setcofactors)


###################################################
### code chunk number 35: ManualAutoStart (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## mqmplot.cofactors(maug, autocofactors, justdots=TRUE)
## mqmplot.cofactors(maug, setcofactors, justdots=TRUE)


###################################################
### code chunk number 36: MQM-tour.Rnw:642-644
###################################################
# plot result of cofactor selection
par(mfrow = c(2,1))
mqmplot.cofactors(maug, autocofactors, justdots=TRUE)
mqmplot.cofactors(maug, setcofactors, justdots=TRUE)


###################################################
### code chunk number 37: ManualAuto (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_backw))
## plot(mqmgetmodel(mqm_auto))


###################################################
### code chunk number 38: MQM-tour.Rnw:660-662
###################################################
# plot result of cofactor backward elimination
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_backw))
plot(mqmgetmodel(mqm_auto))


###################################################
### code chunk number 39: Backward1multi (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_backw))
## plot(mqm_backw)


###################################################
### code chunk number 40: MQM-tour.Rnw:677-679
###################################################
# plot result of cofactor backward elimination
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_backw))
plot(mqm_backw)


###################################################
### code chunk number 41: Backward2 (eval = FALSE)
###################################################
## plot(m_one, mqm_backw, col=c("black","green"), lty=1:2)
## legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)


###################################################
### code chunk number 42: Backward2multi (eval = FALSE)
###################################################
## plot(m_one, mqm_backw, col=c("black","green"), lty=1:2)
## legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)


###################################################
### code chunk number 43: MQM-tour.Rnw:705-706
###################################################
plot(m_one, mqm_backw, col=c("black","green"), lty=1:2)
legend("topleft", c("scanone","MQM"), col=c("black","green"), lwd=1)


###################################################
### code chunk number 44: FigLowAlpha (eval = FALSE)
###################################################
## mqm_backw_low <- mqmscan(maug, setcofactors, cofactor.significance=0.002)
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_backw_low))
## plot(mqm_backw,mqm_backw_low, col=c("blue","green"), lty=1:2)
## legend("topleft", c("Significance=0.02","Significance=0.002"), 
##        col=c("blue","green"), lwd=1)


###################################################
### code chunk number 45: MQM-tour.Rnw:745-746
###################################################
mqm_backw_low <- mqmscan(maug, setcofactors, cofactor.significance=0.002)
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_backw_low))
plot(mqm_backw,mqm_backw_low, col=c("blue","green"), lty=1:2)
legend("topleft", c("Significance=0.02","Significance=0.002"), 
       col=c("blue","green"), lwd=1)


###################################################
### code chunk number 46: AutoCofactor (eval = FALSE)
###################################################
## mqmplot.singletrait(mqm_backw_low, extended=TRUE)


###################################################
### code chunk number 47: MQM-tour.Rnw:775-776
###################################################
mqmplot.singletrait(mqm_backw_low, extended=TRUE)


###################################################
### code chunk number 48: QTLeffects (eval = FALSE)
###################################################
## dirresults <- mqmplot.directedqtl(multitrait, mqm_backw_low)


###################################################
### code chunk number 49: MQM-tour.Rnw:809-810
###################################################
dirresults <- mqmplot.directedqtl(multitrait, mqm_backw_low)


###################################################
### code chunk number 50: MainEffectsD1 (eval = FALSE)
###################################################
## plotPXG(multitrait, marker="GH.117C")


###################################################
### code chunk number 51: MQM-tour.Rnw:828-829
###################################################
plotPXG(multitrait, marker="GH.117C")


###################################################
### code chunk number 52: epistatic1 (eval = FALSE)
###################################################
## effectplot(multitrait, mname1="GH.117C", mname2="GA1")


###################################################
### code chunk number 53: MQM-tour.Rnw:850-851
###################################################
effectplot(multitrait, mname1="GH.117C", mname2="GA1")


###################################################
### code chunk number 54: epistatic2 (eval = FALSE)
###################################################
## effectplot(multitrait, mname1="PVV4", mname2="GH.117C")


###################################################
### code chunk number 55: MQM-tour.Rnw:877-878
###################################################
effectplot(multitrait, mname1="PVV4", mname2="GH.117C")


###################################################
### code chunk number 56: MQM-tour.Rnw:921-924
###################################################
require(snow)
results <- mqmpermutation(maug, scanfunction=mqmscan, cofactors=setcofactors,
                          n.cluster=2, n.perm=25, batchsize=25)


###################################################
### code chunk number 57: MQM-tour.Rnw:927-928
###################################################
mqmplot.permutations(results)


###################################################
### code chunk number 58: MQM-tour.Rnw:939-941
###################################################
resultsrqtl <- mqmprocesspermutation(results)
summary(resultsrqtl)


###################################################
### code chunk number 59: MQM-tour.Rnw:962-965
###################################################
data(multitrait)
m_imp <- fill.geno(multitrait)
mqmscanfdr(m_imp, mqmscanall, cofactors=setcofactors, n.cluster=2)


###################################################
### code chunk number 60: MQM-tour.Rnw:1009-1012
###################################################
data(multitrait)
m_imp <- fill.geno(multitrait)
mqm_imp5 <- mqmscan(m_imp, pheno.col=1:5, n.cluster=2)


###################################################
### code chunk number 61: MQM-tour.Rnw:1015-1016
###################################################
mqmplot.multitrait(mqm_imp5, type="image")


###################################################
### code chunk number 62: MQM-tour.Rnw:1024-1027
###################################################
cofactorlist <- mqmsetcofactors(m_imp, 3)
mqm_imp5 <- mqmscan(m_imp, pheno.col=1:5 , cofactors=cofactorlist,
                    n.cluster=2)


###################################################
### code chunk number 63: MQM-tour.Rnw:1030-1031
###################################################
mqmplot.multitrait(mqm_imp5, type="image")


###################################################
### code chunk number 64: MQM-tour.Rnw:1041-1042
###################################################
mqmplot.multitrait(mqm_imp5, type="lines")


###################################################
### code chunk number 65: MQM-tour.Rnw:1062-1063
###################################################
mqmplot.circle(m_imp, mqm_imp5)


###################################################
### code chunk number 66: MQM-tour.Rnw:1072-1073
###################################################
mqmplot.circle(m_imp, mqm_imp5, highlight=2)


###################################################
### code chunk number 67: MQM-tour.Rnw:1095-1098
###################################################
data(locations)
multiloc <- addloctocross(m_imp, locations)
mqmplot.cistrans(mqm_imp5, multiloc, 5, FALSE, TRUE)


###################################################
### code chunk number 68: MQM-tour.Rnw:1111-1112
###################################################
mqmplot.circle(multiloc, mqm_imp5, highlight=2)


