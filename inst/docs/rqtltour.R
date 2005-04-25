
library(qtl)

ls()

help(read.cross)
?read.cross

data(hyper)
ls()
?hyper

summary(hyper)

nind(hyper)
nphe(hyper)
nchr(hyper)
totmar(hyper)
nmar(hyper)

plot(hyper)

plot.missing(hyper)
plot.map(hyper)
hist(hyper$pheno[,1], breaks=30)

plot.missing(hyper,reorder=TRUE)

hyper <- drop.nullmarkers(hyper)
totmar(hyper)

hyper <- est.rf(hyper)
plot.rf(hyper)
plot.rf(hyper,c(1,4))

plot.rf(hyper,6)
plot.missing(hyper,6)

newmap <- est.map(hyper, error.prob=0.01, verbose=TRUE)
plot.map(hyper, newmap)

hyper <- replace.map(hyper, newmap)

hyper <- calc.genoprob(hyper, error.prob=0.01)
hyper <- calc.errorlod(hyper, error.prob=0.01)

plot.errorlod(hyper)
top.errorlod(hyper)
plot.errorlod(hyper, chr=c(4,11,16))

plot.geno(hyper, chr=16, ind=71:90, min.sep=4)

hyper <- calc.genoprob(hyper, step=2, err=0.01)
plot.info(hyper)
plot.info(hyper, chr=c(1,4,15))
plot.info(hyper, chr=c(1,4,15), method="entropy")
plot.info(hyper, chr=c(1,4,15), method="variance")

out.em <- scanone(hyper)
out.hk <- scanone(hyper, method="hk")

hyper <- sim.geno(hyper, step=2, n.draws=16)
out.imp <- scanone(hyper, method="imp")

summary(out.em)
summary(out.em, 3)
summary(out.hk, 3)
summary(out.imp, 3)

max(out.em)
max(out.hk)
max(out.imp)

plot(out.em, chr=c(1,4,15))
plot(out.hk, out.imp, out.em, chr=c(1,4,15), col=c("red","blue","black"), lty=1)
plot(out.em, chr=c(1,4,15))
plot(out.hk, chr=c(1,4,15), col="blue", add=TRUE)
plot(out.imp, chr=c(1,4,15), col="red", add=TRUE)

operm.hk <- scanone(hyper, method="hk", n.perm=10)
quantile(operm.hk, 0.95)

save.image()

hyper.coarse <- calc.genoprob(hyper, step=10, error.prob=0.01)

out2.hk <- scantwo(hyper.coarse, method="hk")

summary(out2.hk, c(8,3,3))
summary(out2.hk, c(0,4,1000))
summary(out2.hk, c(0,1000,4))

plot(out2.hk)
plot(out2.hk,chr=c(1,4))

max(out2.hk)

chr <- c(1, 1, 4, 6, 15)
pos <- c(50, 76, 30, 70, 20)
qtl <- makeqtl(hyper, chr, pos)

my.formula <- y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q4:Q5
out.fitqtl <- fitqtl(hyper$pheno[,1], qtl, formula=my.formula)
summary(out.fitqtl)

ls()

data(badorder)
summary(badorder)
plot(badorder)

badorder <- est.rf(badorder)
plot.rf(badorder)

plot.rf(badorder, chr=1)

newmap <- est.map(badorder, verbose=TRUE)
plot.map(badorder, newmap)

plot.rf(badorder, chr=2:3)

badorder$geno[[2]]$map[6]
badorder$geno[[3]]$map[5]

badorder <- movemarker(badorder, "D2M937", 3, 48)
badorder <- movemarker(badorder, "D3M160", 2, 28.8)

plot.rf(badorder, chr=2:3)

rip1 <- ripple(badorder, chr=1, window=6)
summary(rip1)

rip2 <- ripple(badorder, chr=1, window=3, err=0.01, method="likelihood")
summary(rip2)

badorder.rev <- switch.order(badorder, 1, rip1[2,])
rip1r <- ripple(badorder.rev, chr=1, window=6)
summary(rip1r)

badorder.rev <- switch.order(badorder.rev, 1, rip1r[2,])
rip2r <- ripple(badorder.rev, chr=1, window=3, err=0.01)
summary(rip2r)

badorder.rev <- est.rf(badorder.rev)
plot.rf(badorder.rev, 1)

data(listeria)
summary(listeria)
plot(listeria)
plot.missing(listeria)

y <- log(listeria$pheno[,1])
listeria$pheno <- cbind(listeria$pheno, logSurv=y)
plot(listeria)

listeria <- est.rf(listeria)
plot.rf(listeria)
plot.rf(listeria,c(5,13))

newmap <- est.map(listeria)
plot.map(listeria, newmap)

listeria <- calc.genoprob(listeria,error.prob=0.01)
listeria <- calc.errorlod(listeria,error.prob=0.01)
plot.errorlod(listeria)
top.errorlod(listeria)
plot.errorlod(listeria,c(5,13))
plot.geno(listeria, chr=13, ind=61:70, min.sep=2 )

listeria <- calc.genoprob(listeria, step=2)
out.2p <- scanone(listeria, pheno.col=2, model="2part", upper=TRUE)

summary(out.2p)
summary(out.2p, 4.5)
plot(out.2p)
plot(out.2p, out.2p[,-3], out.2p[,-(3:4)], chr=c(1,5,13,15),
     lty=1, col=c("black", "red", "blue"))

operm.2p <- scanone(listeria, model="2part", pheno.col=2,
                    upper=TRUE, n.perm=3)
apply(operm.2p, 2, quantile, 0.95)

z <- y <- x <- listeria$pheno[,2]
mx <- max(x, na.rm=TRUE)
y[!is.na(x) & x==mx] <- NA
z[!is.na(x) & x<mx] <- 0
z[!is.na(x) & x==mx] <- 1
listeria$pheno <- cbind(listeria$pheno, logSurv2=y, binary=z)
plot(listeria)

out.mu <- scanone(listeria, pheno.col=3)
plot(out.mu, out.2p[,-(3:4)], chr=c(1,5,13,15))

out.p <- scanone(listeria, pheno.col=4, model="binary")
plot(out.p, out.2p[,-3], chr=c(1,5,13,15))

out.np1 <- scanone(listeria, model="np", ties.random=TRUE)
out.np2 <- scanone(listeria, model="np", ties.random=FALSE)
plot(out.np1, out.np2)
plot(out.2p, out.np1, out.np2, chr=c(1,5,13,15), lty=1,
     col=c("black", "blue", "red"))

data(fake.bc)
summary(fake.bc)
plot(fake.bc)

fake.bc <- calc.genoprob(fake.bc, step=2.5)
out1.nocovar <- scanone(fake.bc, pheno.col=1)
out2.nocovar <- scanone(fake.bc, pheno.col=2)

ac <- fake.bc$pheno[,"age"]
out1.covar.a <- scanone(fake.bc, pheno.col=1, addcov=ac)
out2.covar.a <- scanone(fake.bc, pheno.col=2, addcov=ac)

ac <- fake.bc$pheno[,c("sex","age")]
ic <- fake.bc$pheno[,"sex"]
out1.covar.b <- scanone(fake.bc, pheno.col=1, addcov=ac, intcov=ic)
out2.covar.b <- scanone(fake.bc, pheno.col=2, addcov=ac, intcov=ic)

summary(out1.nocovar, 3)
summary(out1.covar.a, 3)
summary(out1.covar.b, 3)

summary(out2.nocovar, 3)
summary(out2.covar.a, 3)
summary(out2.covar.b, 3)

plot(out1.nocovar, out1.covar.a, out1.covar.b, lty=1,
     chr=c(2,5,10), col=c("black","blue","red"))
plot(out2.nocovar, out2.covar.a, out2.covar.b, lty=1,
     chr=c(2,5,10), col=c("black","blue","red"))

data(fake.bc)

class(fake.bc)

names(fake.bc)

fake.bc$pheno

names(fake.bc$geno)
sapply(fake.bc$geno, class)

names(fake.bc$geno[[3]])
fake.bc$geno[[3]]$data[1:5,]
fake.bc$geno[[3]]$map

names(fake.bc$geno[[3]])
fake.bc <- calc.genoprob(fake.bc, step=10, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- sim.geno(fake.bc, step=10, n.draws=8, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- argmax.geno(fake.bc, step=10, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- calc.errorlod(fake.bc, err=0.01)
names(fake.bc$geno[[3]])

names(fake.bc)
fake.bc <- est.rf(fake.bc)
names(fake.bc)

