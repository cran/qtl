#####################################################################
#
# xchr.R
#
# copyright (c) 2005, Karl W Broman, Johns Hopkins University
# last modified Apr, 2005
# first written Apr, 2004
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: Utilities for dealing with the X chromosome.
#           getsex, getgenonames, reviseXdata, scanoneXnull
#           revisecovar
#           [See also fixXgeno.bc & fixXgeno.f2 in read.cross.R]
#
######################################################################

# get sex and pgm columns from phenotype data
getsex <-
function(cross)
{
  phe.names <- names(cross$pheno)

  sex.column <- grep("^[Ss][Ee][Xx]$", phe.names)
  pgm.column <- grep("^[Pp][Gg][Mm]$", phe.names)

  if(length(sex.column)==0) { # no sex included
    sex <- NULL
  }
  else {
    if(length(sex.column)>1)
      warning("'sex' included multiple times.  Using the first one.")
    temp <- cross$pheno[,sex.column[1]]
    if(is.numeric(temp)) {
      if(any(!is.na(temp) & temp != 0 & temp != 1)) {
        warning("Sex column should be coded as 0=female 1=male; sex ignored.")
        sex <- NULL
      }
      else sex <- temp
    }
    else {
      if(!is.factor(temp)) temp <- as.factor(temp)

      if(length(levels(temp)) == 1) {
        if(levels(temp) == "F" || levels(temp)=="f") sex <- rep(0,nind(cross))
        else if(levels(temp) == "M" || levels(temp)=="m") sex <- rep(1,nind(cross))
        else 
          warning("Sex column should be coded as 0=female 1=male; sex ignored.")
      }
      else if(length(levels(temp)) > 2) {
        warning("Sex column should be coded as a two-level factor; sex ignored.")
        sex <- NULL
      }
      else { # is a factor with two levels
        lev <- levels(temp)
        if(length(grep("^[Ff]",lev))>0 &&
           length(males <- grep("^[Mm]",lev))>0) {
          temp <- as.character(temp)
          sex <- rep(0,length(temp))
          sex[is.na(temp)] <- NA
          sex[!is.na(temp) & temp==lev[males]] <- 1
        }
        else 
          warning("Don't understand levels in sex column; sex ignored.")
      }
    }
  }

  if(length(pgm.column)==0) { # no pgm included
    pgm <- NULL
  }
  else {
    if(length(pgm.column)>1)
      warning("'pgm' included multiple times.  Using the first one.")
    temp <- cross$pheno[,pgm.column[1]]
    if(!is.numeric(temp))
      temp <- as.numeric(temp)-1
    if(any(!is.na(temp) & temp != 0 & temp != 1)) {
      warning("pgm column should be coded as 0/1; pgm ignored.")
      pgm <- NULL
    }
    else pgm <- temp
  }

  list(sex=sex,pgm=pgm)
}
          


# get names of genotypes
# used in discan, effectplot, plot.pxg, scanone, scantwo, vbscan, reviseXdata
getgenonames <-
function(type=c("f2","bc","f2ss","riself","risib","4way"),
         chrtype=c("A","X"), expandX=c("simple","standard","full"),
         sexpgm)
{  
  sex <- sexpgm$sex
  pgm <- sexpgm$pgm

  # get rid of missing sex and pgm values, if there are any
  if(length(sex)>0) sex <- sex[!is.na(sex)]
  if(length(pgm)>0) pgm <- pgm[!is.na(pgm)]

  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  expandX <- match.arg(expandX)

  if(type=="riself" || type=="risib") 
    gen.names <- c("AA","BB")

  else if(type == "4way") {
    if(chrtype=="A") gen.names <- c("AC","BC","AD","BD")
    else gen.names <- c("AC","BC","AY","BY")
  }

  else if(type == "bc") {
    if(chrtype=="A") gen.names <- c("AA","AB") # autosome
    else { # X chromosome
 
#                 simple     standard       full      
#   -both sexes   A-/AB/BY   AA/AB/AY/BY    same as std
#   -all females  AA/AB      same           same
#   -all males    AY/BY      same           same

      if(length(sex)==0 || all(sex==0)) # all females
        gen.names <- c("AA","AB")
      else if(all(sex==1)) # all males
        gen.names <- c("AY","BY")
      else { # some of each
        if(expandX == "simple") gen.names <- c("A-", "AB", "BY")
        else gen.names <- c("AA","AB","AY","BY")
      }
    }
  }

  else { # intercross
    if(chrtype == "A")  # autosomal
      gen.names <- c("AA","AB","BB")
    else { # X chromsome

# both crosses     simple     standard         full
#   -both sexes   A-/AB/B-    AA/AB/BB/AY/BY   AA/AB1/AB2/BB/AY/BY
#   -all females  AA/AB/BB    same as simple   AA/AB1/AB2/BB
#   -all males    AY/BY       same             same
# forw cross
#   -both sexes   A-/AB/BY    AA/AB/AY/BY      same as std
#   -all females  AA/AB       same             same
#   -all males    AY/BY       same             same
# backw cross
#   -both sexes   B-/AB/AY    BB/AB/AY/BY      same as std
#   -all females  BB/AB       same             same
#   -all males    AY/BY       same             same

      if(length(sex)==0 || all(sex==0)) { # all females
        if(length(pgm)==0 || all(pgm==0)) # all forw dir
          gen.names <- c("AA","AB")
        else if(all(pgm==1))  # all backw dir
          gen.names <- c("BB","AB")
        else { # some of each direction
          if(expandX=="full") gen.names <- c("AA","ABf","ABr","BB")
          else gen.names <- c("AA","AB","BB")
        }
      }
      else if(all(sex==1))  # all males
        gen.names <- c("AY","BY")
      else { # some of each sex
        if(length(pgm)==0 || all(pgm==0)) { # all forw
          if(expandX=="simple") gen.names <- c("A-","AB","BY")
          else gen.names <- c("AA","AB","AY","BY")
        }
        else if (all(pgm==1)) { # all backw
          if(expandX=="simple") gen.names <- c("B-","AB","AY")
          else gen.names <- c("BB","AB","AY","BY")
        }
        else { # some of each dir
          if(expandX=="simple") gen.names <- c("A-","AB","B-")
          else if(expandX=="standard")
            gen.names <- c("AA","AB","BB","AY","BY")
          else
            gen.names <- c("AA","ABf","ABr","BB","AY","BY")
        }
      }
    }
  }

  gen.names
}

# revise genotype data, probabilities or imputations for the X chromosome
reviseXdata <-
function(type=c("f2ss","f2","bc"), expandX=c("simple","standard","full"),
         sexpgm, geno, prob, draws, pairprob)
{
  type <- match.arg(type)
  expandX <- match.arg(expandX)

  sex <- sexpgm$sex
  pgm <- sexpgm$pgm

  notmissing <- (!missing(geno)) + (!missing(prob)) + (!missing(draws)) +
      (!missing(pairprob))
  if(notmissing == 0)
    stop("Provide one of geno, prob, draws, pairprob.")
  if(notmissing > 1)
    stop("Provide just one of geno, prob, draws, pairprob.")

  # get genonames
  genonames <- getgenonames(type, "X", expandX, sexpgm)

  if(type == "bc") { # backcross

    if(length(sex)==0 || all(sex==0) || all(sex==1)) { # all one sex
      # no changes necessary
      if(!missing(geno)) return(geno)
      else if(!missing(prob)) {
        dimnames(prob)[[3]] <- genonames
        return(prob)
      }
      else if(!missing(draws)) 
        return(draws)
      else # pairprob
        return(pairprob)
    }

    else { # both sexes

      if(!missing(geno)) {
        gmale <- geno[sex==1,]
        if(expandX=="simple") 
          gmale[!is.na(gmale) & gmale==2] <- 3
        else {
          gmale[!is.na(gmale) & gmale==1] <- 3
          gmale[!is.na(gmale) & gmale==2] <- 4
        }
        geno[sex==1,] <- gmale
        return(geno)
      }

      else if(!missing(draws)) {
        gmale <- draws[sex==1,,]
        if(expandX=="simple") 
          gmale[gmale==2] <- 3
        else {
          gmale[gmale==1] <- 3
          gmale[gmale==2] <- 4
        }
        draws[sex==1,,] <- gmale
        return(draws)
      }

      else if(!missing(prob)) {
        dimprob <- dim(prob)
        dimprob[3] <- length(genonames)
        newprob <- array(0,dim=dimprob)
        dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
        newprob[sex==0,,1:2] <- prob[sex==0,,]

        if(expandX=="simple") {
          newprob[sex==1,,1] <- prob[sex==1,,1]
          newprob[sex==1,,3] <- prob[sex==1,,2]
        }
        else {
          newprob[sex==1,,3] <- prob[sex==1,,1]
          newprob[sex==1,,4] <- prob[sex==1,,2]
        }
        return(newprob)
      }

      else { # pairprob
        dimpairprob <- dim(pairprob)
        dimpairprob[3] <- dimpairprob[4] <- length(genonames)
        newpairprob <- array(0,dim=dimpairprob)
        newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]
        
        if(expandX=="simple") {
          newpairprob[sex==1,,1,1] <- pairprob[sex==1,,1,1]
          newpairprob[sex==1,,1,3] <- pairprob[sex==1,,1,2]
          newpairprob[sex==1,,3,1] <- pairprob[sex==1,,2,1]
          newpairprob[sex==1,,3,3] <- pairprob[sex==1,,2,2]
        }
        else {
          newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
          newpairprob[sex==1,,3,4] <- pairprob[sex==1,,1,2]
          newpairprob[sex==1,,4,3] <- pairprob[sex==1,,2,1]
          newpairprob[sex==1,,4,4] <- pairprob[sex==1,,2,2]
        }
        return(newpairprob)
      }
          
    } # end of "both sexes" / backcross

  } # end of backcross

  else { # intercross

    if(length(sex)==0 || all(sex==0)) { # all females

      if(length(pgm)==0 || all(pgm==0) || all(pgm==1)) { # one dir, females
        if(!missing(geno)) return(geno)
        else if(!missing(draws)) return(draws)
        else if(!missing(pairprob)) return(pairprob)
        else {
          dimnames(prob)[[3]] <- genonames
          return(prob)
        }
      }
        
      else { # both dir, females
        if(!missing(geno)) {
          gback <- geno[pgm==1,]
          gback[!is.na(gback) & gback==1] <- 3
          geno[pgm==1,] <- gback
          return(geno)
        }
        else if(!missing(draws)) {
          gback <- draws[pgm==1,,]
          gback[!is.na(gback) & gback==1] <- 3
          draws[pgm==1,,] <- gback
          return(draws)
        }
        else if(!missing(prob)) {
          dimprob <- dim(prob)
          dimprob[3] <- length(genonames)
          newprob <- array(0,dim=dimprob)
          dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
          newprob[pgm==0,,1:2] <- prob[pgm==0,,]

          if(expandX!="full") { # simple/standard
            newprob[pgm==1,,3] <- prob[pgm==1,,1]
            newprob[pgm==1,,2] <- prob[pgm==1,,2]
          }
          else {
            newprob[pgm==1,,4] <- prob[pgm==1,,1]
            newprob[pgm==1,,3] <- prob[pgm==1,,2]
          }
          return(newprob)
        }
        else { # pairprob
          dimpairprob <- dim(pairprob)
          dimpairprob[3] <- dimpairprob[4] <- length(genonames)
          newpairprob <- array(0,dim=dimpairprob)
          newpairprob[pgm==0,,1:2,1:2] <- pairprob[pgm==0,,,]
        
          if(expandX!="full") { # simple/standard
            newpairprob[pgm==1,,3,3] <- pairprob[pgm==1,,1,1]
            newpairprob[pgm==1,,3,2] <- pairprob[pgm==1,,1,2]
            newpairprob[pgm==1,,2,3] <- pairprob[pgm==1,,2,1]
            newpairprob[pgm==1,,2,2] <- pairprob[pgm==1,,2,2]
          }
          else {
            newpairprob[pgm==1,,4,4] <- pairprob[pgm==1,,1,1]
            newpairprob[pgm==1,,4,3] <- pairprob[pgm==1,,1,2]
            newpairprob[pgm==1,,3,4] <- pairprob[pgm==1,,2,1]
            newpairprob[pgm==1,,3,3] <- pairprob[pgm==1,,2,2]
          }
        return(newpairprob)
        }
      }
    }
    else if(all(sex==1))  { # all males
      if(!missing(geno)) return(geno)
      else if(!missing(draws)) return(draws)
      else if(!missing(pairprob)) return(pairprob)
      else {
        dimnames(prob)[[3]] <- genonames
        return(prob)
      }
    }

    else { # both sexes

      if(length(pgm)==0 || all(pgm==0)) { # both sexes, forw dir
        if(!missing(geno)) {
          gmale <- geno[sex==1,]
          if(expandX!="full") 
            gmale[!is.na(gmale) & gmale==2] <- 3
          else {
            gmale[!is.na(gmale) & gmale==1] <- 3
            gmale[!is.na(gmale) & gmale==2] <- 4
          }
          geno[sex==1,] <- gmale
          return(geno)
        }

        else if(!missing(draws)) {
          gmale <- draws[sex==1,,]
          if(expandX!="full") 
            gmale[gmale==2] <- 3
          else {
            gmale[gmale==1] <- 3
            gmale[gmale==2] <- 4
          }
          draws[sex==1,,] <- gmale
          return(draws)
        }

        else if(!missing(prob)) {
          dimprob <- dim(prob)
          dimprob[3] <- length(genonames)
          newprob <- array(0,dim=dimprob)
          dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
          newprob[sex==0,,1:2] <- prob[sex==0,,]

          if(expandX=="simple") {
            newprob[sex==1,,1] <- prob[sex==1,,1]
            newprob[sex==1,,3] <- prob[sex==1,,2]
          }
          else {
            newprob[sex==1,,3] <- prob[sex==1,,1]
            newprob[sex==1,,4] <- prob[sex==1,,2]
          }
          return(newprob)
        }

        else { # pairprob
          dimpairprob <- dim(pairprob)
          dimpairprob[3] <- dimpairprob[4] <- length(genonames)
          newpairprob <- array(0,dim=dimpairprob)
          newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]
        
          if(expandX=="simple") {
            newpairprob[sex==1,,1,1] <- pairprob[sex==1,,1,1]
            newpairprob[sex==1,,1,3] <- pairprob[sex==1,,1,2]
            newpairprob[sex==1,,3,1] <- pairprob[sex==1,,2,1]
            newpairprob[sex==1,,3,3] <- pairprob[sex==1,,2,2]
          }
          else {
            newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
            newpairprob[sex==1,,3,4] <- pairprob[sex==1,,1,2]
            newpairprob[sex==1,,4,3] <- pairprob[sex==1,,2,1]
            newpairprob[sex==1,,4,4] <- pairprob[sex==1,,2,2]
          }
          return(newpairprob)
        }
      } # both sexes, forw dir

      if(all(pgm==1)) { # both sexes, backw dir
        if(!missing(geno)) {
          gmale <- geno[sex==1,]
          if(expandX!="full") {
            gmale[!is.na(gmale) & gmale==1] <- 3
            gmale[!is.na(gmale) & gmale==2] <- 1
          }
          else {
            gmale[!is.na(gmale) & gmale==1] <- 3
            gmale[!is.na(gmale) & gmale==2] <- 4
          }
          geno[sex==1,] <- gmale
          return(geno)
        }

        else if(!missing(draws)) {
          gmale <- draws[sex==1,,]
          if(expandX!="full") {
            gmale[gmale==1] <- 3
            gmale[gmale==2] <- 1
          }
          else {
            gmale[gmale==1] <- 3
            gmale[gmale==2] <- 4
          }
          draws[sex==1,,] <- gmale
          return(draws)
        }

        else if(!missing(prob)) {
          dimprob <- dim(prob)
          dimprob[3] <- length(genonames)
          newprob <- array(0,dim=dimprob)
          dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
          newprob[sex==0,,1:2] <- prob[sex==0,,]

          if(expandX=="simple") {
            newprob[sex==1,,3] <- prob[sex==1,,1]
            newprob[sex==1,,1] <- prob[sex==1,,2]
          }
          else {
            newprob[sex==1,,3] <- prob[sex==1,,1]
            newprob[sex==1,,4] <- prob[sex==1,,2]
          }
          return(newprob)
        }

        else { # pairprob
          dimpairprob <- dim(pairprob)
          dimpairprob[3] <- dimpairprob[4] <- length(genonames)
          newpairprob <- array(0,dim=dimpairprob)
          newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]
        
          if(expandX=="simple") {
            newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
            newpairprob[sex==1,,1,3] <- pairprob[sex==1,,2,1]
            newpairprob[sex==1,,3,1] <- pairprob[sex==1,,1,2]
            newpairprob[sex==1,,1,1] <- pairprob[sex==1,,2,2]
          }
          else {
            newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
            newpairprob[sex==1,,3,4] <- pairprob[sex==1,,1,2]
            newpairprob[sex==1,,4,3] <- pairprob[sex==1,,2,1]
            newpairprob[sex==1,,4,4] <- pairprob[sex==1,,2,2]
          }
          return(newpairprob)
        }
      } # both sexes, backw dir

      else { # both dir, both sexes

        if(!missing(geno)) {
          gmale <- geno[sex==1,]
          gfemaler <- geno[sex==0 & pgm==1,]
          if(expandX=="simple") {
            gmale[!is.na(gmale) & gmale==2] <- 3
            gfemaler[!is.na(gfemaler) & gfemaler==1] <- 3
          }
          else if(expandX=="standard") {
            gmale[!is.na(gmale) & gmale==1] <- 4
            gmale[!is.na(gmale) & gmale==2] <- 5
            gfemaler[!is.na(gfemaler) & gfemaler==1] <- 3
          }
          else {
            gmale[!is.na(gmale) & gmale==1] <- 5
            gmale[!is.na(gmale) & gmale==2] <- 6
            gfemaler[!is.na(gfemaler) & gfemaler==1] <- 4
            gfemaler[!is.na(gfemaler) & gfemaler==2] <- 3
          }
          geno[sex==1,] <- gmale
          geno[sex==0 & pgm==1,] <- gfemaler
          return(geno)
        }

        else if(!missing(draws)) {
          gmale <- draws[sex==1,,]
          gfemaler <- draws[sex==0 & pgm==1,,]
          if(expandX=="simple") {
            gmale[gmale==2] <- 3
            gfemaler[gfemaler==1] <- 3
          }
          else if(expandX=="standard") {
            gmale[gmale==1] <- 4
            gmale[gmale==2] <- 5
            gfemaler[gfemaler==1] <- 3
          }
          else {
            gmale[gmale==1] <- 5
            gmale[gmale==2] <- 6
            gfemaler[gfemaler==1] <- 4
            gfemaler[gfemaler==2] <- 3
          }
          draws[sex==1,,] <- gmale
          draws[sex==0 & pgm==1,,] <- gfemaler
          return(draws)
        }

        else if(!missing(prob)) {
          dimprob <- dim(prob)
          dimprob[3] <- length(genonames)
          newprob <- array(0,dim=dimprob)
          dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
          newprob[sex==0 & pgm==0,,1:2] <- prob[sex==0 & pgm==0,,]

          if(expandX=="simple") {
            newprob[sex==1,,1] <- prob[sex==1,,1]
            newprob[sex==1,,3] <- prob[sex==1,,2]
            newprob[sex==0 & pgm==1,,3] <- prob[sex==0 & pgm==1,,1]
            newprob[sex==0 & pgm==1,,2] <- prob[sex==0 & pgm==1,,2]
          }
          else if(expandX=="standard") {
            newprob[sex==1,,4] <- prob[sex==1,,1]
            newprob[sex==1,,5] <- prob[sex==1,,2]
            newprob[sex==0 & pgm==1,,3] <- prob[sex==0 & pgm==1,,1]
            newprob[sex==0 & pgm==1,,2] <- prob[sex==0 & pgm==1,,2]
          }
          else {
            newprob[sex==1,,5] <- prob[sex==1,,1]
            newprob[sex==1,,6] <- prob[sex==1,,2]
            newprob[sex==0 & pgm==1,,4] <- prob[sex==0 & pgm==1,,1]
            newprob[sex==0 & pgm==1,,3] <- prob[sex==0 & pgm==1,,2]
          }
          return(newprob)
        }

        else { # pairprob
          dimpairprob <- dim(pairprob)
          dimpairprob[3] <- dimpairprob[4] <- length(genonames)
          newpairprob <- array(0,dim=dimpairprob)
          newpairprob[sex==0 & pgm==0,,1:2,1:2] <- pairprob[sex==0 & pgm==0,,,]
        
          male <- (sex==1)
          femaler <- (sex==0) & (pgm==1)
          if(expandX=="simple") {
            newpairprob[male,,1,1] <- pairprob[male,,1,1]
            newpairprob[male,,1,3] <- pairprob[male,,1,2]
            newpairprob[male,,3,1] <- pairprob[male,,2,1]
            newpairprob[male,,3,3] <- pairprob[male,,2,2]

            newpairprob[femaler,,3,3] <- pairprob[femaler,,1,1]
            newpairprob[femaler,,3,2] <- pairprob[femaler,,1,2]
            newpairprob[femaler,,2,3] <- pairprob[femaler,,2,1]
            newpairprob[femaler,,2,2] <- pairprob[femaler,,2,2]
          }
          else if(expandX=="standard") {
            newpairprob[male,,4,4] <- pairprob[male,,1,1]
            newpairprob[male,,4,5] <- pairprob[male,,1,2]
            newpairprob[male,,5,4] <- pairprob[male,,2,1]
            newpairprob[male,,5,5] <- pairprob[male,,2,2]

            newpairprob[femaler,,3,3] <- pairprob[femaler,,1,1]
            newpairprob[femaler,,3,2] <- pairprob[femaler,,1,2]
            newpairprob[femaler,,2,3] <- pairprob[femaler,,2,1]
            newpairprob[femaler,,2,2] <- pairprob[femaler,,2,2]
          }
          else {
            newpairprob[male,,5,5] <- pairprob[male,,1,1]
            newpairprob[male,,5,6] <- pairprob[male,,1,2]
            newpairprob[male,,6,5] <- pairprob[male,,2,1]
            newpairprob[male,,6,6] <- pairprob[male,,2,2]

            newpairprob[femaler,,4,4] <- pairprob[femaler,,1,1]
            newpairprob[femaler,,4,3] <- pairprob[femaler,,1,2]
            newpairprob[femaler,,3,4] <- pairprob[femaler,,2,1]
            newpairprob[femaler,,3,3] <- pairprob[femaler,,2,2]
          }
          return(newpairprob)
        }

      } 
    } 

  } # end of intercross

}

######################################################################
# scanoneXnull
#
# figure out null hypothesis business for scanone on X chromosome
######################################################################
scanoneXnull <-
function(type, sexpgm)
{
  sex <- sexpgm$sex
  pgm <- sexpgm$pgm

  if(type == "f2ss") type <- "f2"
  if(type == "risib" || type=="riself") type <- "bc"

  ### first figure out sex/pgm pattern

  # sex
  if(length(sex)==0 || all(sex==0)) { # all female
    onesex <- allfemale <- TRUE
  }
  else if(all(sex==1)) { # all male
    onesex <- TRUE
    allfemale <- FALSE
  }
  else { # both sexes
    onesex <- allfemale <- FALSE
  }
  # pgm
  if(length(pgm)==0 || all(pgm==0) || all(pgm==1)) # one direction
    onedir <- TRUE
  else onedir <- FALSE

  allmale <- onesex && !allfemale
  bothsex <- !onesex
  bothdir <- !onedir


  ### now figure out the null hypothesis and pull out appropriate
  ### covariates for the null

  # backcross, one sex
  # OR intercross, one dir and one sex
  # OR intercross, both dir and all male
  if((type=="bc" && onesex) ||
     (type=="f2" && ((onedir && onesex) || (bothdir && allmale)))) {
    adjustX <- FALSE
    dfX <- 1
    sexpgmcovar <- sexpgmcovar.alt <- NULL
  }

  # backcross, both sexes
  # OR intercross, one direction and both sexes
  else if((type=="bc" && bothsex) ||
          (type=="f2" && onedir && bothsex)) {
    adjustX <- TRUE
    dfX <- 2
    sexpgmcovar <- cbind(sex)
    sexpgmcovar.alt <- sex+1
  }

  # intercross, both dir and all female
  else if(type=="f2" && bothdir && allfemale) {
    adjustX <- TRUE
    dfX <- 2
    sexpgmcovar <- cbind(pgm)
    sexpgmcovar.alt <- pgm+1
  }

  # intercross, both dir and both sexes
  else {
    adjustX <- TRUE
    dfX <- 3
    sexpgmcovar <- cbind(sex,as.numeric(sex==0 & pgm==1))
    sexpgmcovar.alt <- rep(3,length(sex))
    sexpgmcovar.alt[sex==0 & pgm==0] <- 1
    sexpgmcovar.alt[sex==0 & pgm==1] <- 2
  }

  list(adjustX=adjustX, dfX=dfX, sexpgmcovar=sexpgmcovar,
       sexpgmcovar.alt=sexpgmcovar.alt)
}

######################################################################
# revisecovar
#
# Drop sex and pgm and their interxn as covariates for the X chr.
######################################################################
revisecovar <-
function(sexpgm, covar)
{
  if(is.null(covar) || (is.null(sexpgm$sex) && is.null(sexpgm$pgm)))
    return(covar)

  covar <- as.matrix(covar)

  X <- cbind(1,sexpgm$sex,sexpgm$pgm, sexpgm$sex*sexpgm$pgm)
  nc <- ncol(X)
  keep <- rep(TRUE,ncol(covar))
  for(i in 1:ncol(covar)) {
    if(qr(cbind(X,covar[,i]))$rank <= nc)
      keep[i] <- FALSE
  }
  if(!any(keep)) return(NULL)
  covar[,keep,drop=FALSE]
}

# end of xchr.R
