######################################################################
#
# read.cross.R
#
# copyright (c) 2000-2001, Karl W Broman, Johns Hopkins University
# Oct, 2001; Sept, 2001; Feb, 2001; Sept, 2000; Aug, 2000 
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/qtl package
# Contains: read.cross, read.cross.karl, read.cross.mm,
#           read.cross.gary, read.cross.csv
#
######################################################################


######################################################################
#
# read.cross: read data from an experimental cross
#
######################################################################

read.cross <-
function(format=c("csv","mm","gary","karl"),...)
{
  format <- match.arg(format)
  
  if(format=="karl") # karl's format
    cross <- read.cross.karl(...)
  if(format=="mm") # mapmaker format
    cross <- read.cross.mm(...)
  if(format=="gary") # gary's format
    cross <- read.cross.gary(...)
  if(format=="csv") # comma-delimited format
    cross <- read.cross.csv(...)
  
  cross
}




######################################################################
#
# read.cross.karl
#
# read data in Karl's format
#
######################################################################

read.cross.karl <-
function(dir,genfile,mapfile,phefile)  
{
  # create file names 
  if(missing(genfile)) genfile <- "gen.txt"
  if(missing(mapfile)) mapfile <- "map.txt"
  if(missing(phefile)) phefile <- "phe.txt"

  if(!missing(dir)) {
    # remove ending "/" if it exists
    n <- nchar(dir)
    if(substr(dir,n,n) == "/")
      dir <- substr(dir,0,n-1)
    genfile <- paste(dir,genfile, sep="/")
    mapfile <- paste(dir,mapfile, sep="/")
    phefile <- paste(dir,phefile, sep="/")
  }

  # read data
  geno <- as.matrix(read.table(genfile,na.strings="0"))
  pheno <- as.matrix(read.table(phefile,na.strings="-",header=TRUE))
  tempmap <- scan(mapfile, what=character(),quiet=TRUE)

  # fix up map information
  # number of chromosomes
  n.chr <- as.numeric(tempmap[1])
  n.mar <- 1:n.chr
  g <- map <- geno.data <- vector("list", n.chr)
  cur <- 2
  min.mar <- 1
  names(g) <- as.character(1:n.chr)
  for(i in 1:n.chr) { # loop over chromosomes
    # number of markers
    n.mar[i] <- as.numeric(tempmap[cur])
    cur <- cur+1

    # pull out appropriate portion of genotype data
    geno.data[[i]] <- geno[,min.mar:(min.mar+n.mar[i]-1)]
    min.mar <- min.mar + n.mar[i]

    # recombination fractions
    r <- as.numeric(tempmap[cur:(cur+n.mar[i]-2)])

    # convert to cM distances (w/ Kosambi map function)
    d <- 0.25*log((1+2*r)/(1-2*r))*100

    # convert to locations
    map[[i]] <- round(c(0,cumsum(d)),2)
    cur <- cur+n.mar[i]-1

    # marker names
    names(map[[i]]) <- tempmap[cur:(cur+n.mar[i]-1)]
    dimnames(geno.data[[i]]) <- list(NULL, names(map[[i]]))
    cur <- cur+n.mar[i]

    g[[i]] <- list(data=geno.data[[i]],map=map[[i]])

    # attempt to pull out chromosome number
    mar.names <- names(map[[i]])
    twodig <- grep("[Dd][1-9][0-9][Mm]", mar.names)
    onedig <- grep("[Dd][1-9][Mm]", mar.names)
    xchr <- grep("[Dd][Xx][Mm]", mar.names)

    chr.num <- NULL
    if(length(twodig) > 0)
      chr.num <- c(chr.num,substr(mar.names[twodig],2,3))
    if(length(onedig) > 0)
      chr.num <- c(chr.num,substr(mar.names[onedig],2,2))
    if(length(xchr) > 0)
      chr.num <- c(chr.num,rep("X",length(xchr)))

    # no marker names of the form above
    if(is.null(chr.num)) {
      chr.num <- length(mar.names)
      names(chr.num) <- "1"
    }
    else {
      chr.num <- table(chr.num)
    }
  
    m <- max(chr.num)
    if(m > sum(chr.num)/2 && m > 1) 
      names(g)[i] <- names(chr.num)[chr.num==m][1]

    if(names(g)[i] == "X" || names(g)[i] == "x") class(g[[i]]) <- "X"
    else class(g[[i]]) <- "A"
  }

  # check that data dimensions match 
  n.mar1 <- sapply(g,function(a) ncol(a$data))
  n.mar2 <- sapply(g,function(a) length(a$map))
  n.phe <- ncol(pheno)
  n.ind1 <- nrow(pheno)
  n.ind2 <- sapply(g,function(a) nrow(a$data))
  if(any(n.ind1 != n.ind2)) {
    print(c(n.ind1,n.ind2))
    stop("Number of individuals in genotypes and phenotypes do not match.");
  }
  if(any(n.mar1 != n.mar2)) {
    print(c(n.mar,n.mar2))
    stop("Numbers of markers in genotypes and marker names files do not match.");
  }

  # print some information about the amount of data read
  cat(" --Read the following data:\n");
  cat("\t", n.ind1, " individuals\n");
  cat("\t", sum(n.mar1), " markers\n");
  cat("\t", n.phe, " phenotypes\n");

  # add phenotype names, if missing
  if(is.null(colnames(pheno))) 
    dimnames(pheno) <- list(NULL, paste("phenotype", 1:n.phe,sep=""))

  # determine map type: f2 or bc or 4way?
  if(max(geno[!is.na(geno)])<=2) type <- "bc"
  else if(max(geno[!is.na(geno)])<=5) type <- "f2"
  else type <- "4way"
  cross <- list(geno=g,pheno=pheno)
  class(cross) <- c(type,"cross")

  # check that nothing is strange in the genotype data
  cross.type <- class(cross)[1]
  if(cross.type=="f2") max.gen <- 5
  else if(cross.type=="bc") max.gen <- 2
  else max.gen <- 10

  u <- unique(geno)
  if(any(!is.na(u) & (u > max.gen | u < 1))) 
    stop(paste("There are stange values in the genotype data :",
               paste(u,collapse=":"), "."))

  cross
}
  

######################################################################
#
# read.cross.mm: read data from an experimental cross in mapmaker
#                format.
#
# We need two files: a "raw" file containing the genotype and
# phenotype data and a "map" file containing the chromosomes
# assignments and (optionally) map positions.
#
# The map file contains two or three columns, separated by white
# space, with the chromosome number, marker name (with markers in
# order along the chromosomes) and (optionally) the map position.
#
######################################################################

read.cross.mm <-
function(dir,rawfile,mapfile,estimate.map=TRUE)
{
  # create file names
  if(missing(mapfile))
    stop("You must specify the mapfile containing the chromosome assignments of markers.")
  if(missing(rawfile))
    stop("You must specify the rawfile containing the genotype and phenotype data.")
  if(!missing(dir)) {
    # remove ending "/" if it exists
    n <- nchar(dir)
    if(substring(dir,n,n) == "/")
      dir <- substring(dir,0,n-1)
    mapfile <- paste(dir, mapfile, sep="/")
    rawfile <- paste(dir, rawfile, sep="/")
  }

  # count lines in rawfile
  n.lines <- length(scan(rawfile, what=character(), skip=0, nlines=0,
                         blank.lines.skip=FALSE,quiet=TRUE))
  map <- read.table(mapfile,header=FALSE,as.is=TRUE,blank=FALSE)
  o <- (1:nrow(map))[map[,1]==""]
  if(length(o) > 0) map <- map[-o,]

  # remove any leading *'s from the marker names
  g <- grep("^*",map[,2],extended=FALSE)
  if(length(g) > 0) 
    map[g,2] <- substr(map[g,2],2,nchar(map[g,2]))

  cur.mar <- 0
  cur.phe <- 0

  NEW.symb <- c("1","2","3","4","5","0")
  OLD.symb <- c("A","H","B","D","C","-")

  flag <- 0
  for(i in 1:n.lines) {
    a <- scan(rawfile,what=character(),skip=i-1,nlines=1,
              blank.lines.skip=TRUE,quiet=TRUE)

    if(length(a) == 0) next

    if(length(grep("#", a[1])) != 0) next

    if(flag == 0) {
      flag <- 1
      type <- a[length(a)]
      if(type == "intercross") type <- "f2"
      else if(type == "backcross") type <- "bc"
      else
        stop(paste("File indicates invalid cross type: ", type,
                   ".", sep=""))
    }
    else if(flag == 1) {
      flag <- 2
      n.ind <- as.numeric(a[1])
      n.mar <- as.numeric(a[2])
      n.phe <- as.numeric(a[3])
      cat("\tType of cross:         ", type, "\n")
      cat("\tNumber of individuals: ", n.ind, "\n")
      cat("\tNumber of markers:     ", n.mar, "\n")
      cat("\tNumber of phenotypes:  ", n.phe, "\n")
      convert.symb <- FALSE
      if(length(a) > 3 && a[4] == "symbols") {
        b <- a[-(1:4)]
        convert.symb <- TRUE
        new.symb <- substring(b,1,1)
        old.symb <- substring(b,3,3)
      }

      marnames <- rep("", n.mar)
      geno <- matrix(0,ncol=n.mar,nrow=n.ind)
      if(n.phe == 0) {
        pheno <- matrix(1:n.ind,ncol=1)
        phenames <- c("number")
      }
      else {
        pheno <- matrix(0,ncol=n.phe,nrow=n.ind)
        phenames <- rep("", n.phe)
      }

    }
    else {

      if(substring(a[1],1,1) == "*") {
        cur.mar <- cur.mar+1
        cur.row <- 1

        if(cur.mar > n.mar) { # now reading phenotypes
          cur.phe <- cur.phe+1
          if(cur.phe > n.phe) next 
          phenames[cur.phe] <- substring(a[1],2)
          p <- a[-1]
          p[p=="-"] <- NA
          n <- length(p)
          pheno[cur.row+(0:(n-1)),cur.phe] <- as.numeric(p)
          cur.row <- cur.row + n
        }

        else { # reading genotypes
          marnames[cur.mar] <- substring(a[1],2)
          g <- paste(a[-1],collapse="")
          g <- unlist(strsplit(g,""))

          if(convert.symb) {
            h <- g
            for(j in 1:length(new.symb)) {
              if(any(h==new.symb[j]))
                g[h==new.symb[j]] <- old.symb[j]
            }
          }

          h <- g
          for(j in 1:length(NEW.symb)) {
            if(any(h==OLD.symb[j]))
              g[h==OLD.symb[j]] <- NEW.symb[j]
          }

          n <- length(g)

          geno[cur.row+(0:(n-1)),cur.mar] <- as.numeric(g)
          cur.row <- cur.row + n
        }

      }
      else { # continuation lines
        if(cur.mar > n.mar) { # now reading phenotypes
          a[a=="-"] <- NA
          n <- length(a)
          pheno[cur.row+(0:(n-1)),cur.phe] <- as.numeric(a)
          cur.row <- cur.row + n
        }
        else {
          g <- paste(a,collapse="")
          g <- unlist(strsplit(g,""))
          if(convert.symb) {
            h <- g
            for(j in 1:length(new.symb)) {
              if(any(h==new.symb[j]))
                g[h==new.symb[j]] <- old.symb[j]
            }
          }
          h <- g
          for(j in 1:length(NEW.symb)) {
            if(any(h==OLD.symb[j]))
              g[h==OLD.symb[j]] <- NEW.symb[j]
          }
          n <- length(g)
          geno[cur.row+(0:(n-1)),cur.mar] <- as.numeric(g)
          cur.row <- cur.row + n
        }

      }          

    }
  }
  dimnames(pheno) <- list(NULL, phenames)

  # parse map file
  includes.pos <- FALSE
  if(ncol(map) == 3) includes.pos <- TRUE
  
  chr <- as.character(map[,1])
  markers <- map[,2]
  if(includes.pos) pos <- map[,3]

  Geno <- vector("list",length(unique(chr)))
  names(Geno) <- unique(chr)

  for(i in unique(chr)) {
    mar <- markers[chr == i]

    # create map
    if(includes.pos) map <- pos[chr == i]
    else map <- seq(0,by=5,length=length(mar))
    names(map) <- mar
    
    # pull out genotype data
    o <- match(mar,marnames)
    if(any(is.na(o))) {
      stop(paste("Cannot find markers in genotype data: ",
           paste(mar[is.na(o)],collapse=" "), ".",sep=""))
    }

    if(length(o)==1) data <- matrix(geno[,o],ncol=1)
    else data <- geno[,o]
    # add marker names to data
    colnames(data) <- mar
    # changes 0's to NA's
    data[!is.na(data) & data==0] <- NA

    Geno[[i]] <- list(data=data,map=map)
    if(i=="X" || i=="x") class(Geno[[i]]) <- "X"
    else class(Geno[[i]]) <- "A"
  }

  cross <- list(geno=Geno,pheno=pheno)
  class(cross) <- c(type,"cross")

  if(!includes.pos && estimate.map) {
    cat(" --Estimating genetic map\n")
    newmap <- est.map(cross)
    cross <- replace.map(cross,newmap)
  }

  cross
}


######################################################################
#
# read.cross.gary
#
# read data in Gary's format
#
######################################################################

read.cross.gary <-
function(dir,genfile,mnamesfile,chridfile,phefile,pnamesfile,mapfile)
{
  # create file names 
  if(missing(genfile)) genfile <- "geno.dat"
  if(missing(mnamesfile)) mnamesfile <- "mnames.txt"
  if(missing(chridfile)) chridfile <- "chrid.dat"
  if(missing(phefile)) phefile <- "pheno.dat"
  if(missing(pnamesfile)) pnamesfile <- "pnames.txt"
  if(missing(mapfile)) mapfile <- "markerpos.txt"

  if(!missing(dir)) {
    # remove ending "/" if it exists
    n <- nchar(dir)
    if(substr(dir,n,n) == "/")
      dir <- substr(dir,0,n-1)
    genfile <- paste(dir,genfile, sep="/")
    mnamesfile <- paste(dir,mnamesfile, sep="/")
    chridfile <- paste(dir,chridfile, sep="/")
    phefile <- paste(dir,phefile, sep="/")
    pnamesfile <- paste(dir,pnamesfile, sep="/")
    mapfile <- paste(dir,mapfile, sep="/")
  }

  # read data
  allgeno <- as.matrix(read.table(genfile,na.strings="9"))+1
  pheno <- as.matrix(read.table(phefile,na.strings="-",header=FALSE))
  chr <- scan(chridfile,what=character(),quiet=TRUE)
  mnames <- scan(mnamesfile,what=character(),quiet=TRUE)
  map <- read.table(mapfile,row.names=1)
  pnames <- scan(pnamesfile,what=character(),quiet=TRUE)

  map <- map[mnames,1]

  # fix up map information
  # number of chromosomes
  uchr <- unique(chr)
  n.chr <- length(uchr)
  geno <- vector("list",n.chr)
  names(geno) <- uchr
  min.mar <- 1
  for(i in 1:n.chr) { # loop over chromosomes
    # create map
    temp.map <- map[chr==uchr[i]]

    # deal with any markers that didn't appear in the marker pos file
    if(any(is.na(temp.map))) {
      o <- (1:length(temp.map))[is.na(temp.map)]
      for(j in o) {
        if(j==1 || all(is.na(temp.map[1:(j-1)]))) {
          z <- min((1:length(temp.map))[-o])
          temp.map[j] <- min(temp.map,na.rm=TRUE)-(z-j+1)
        }
        else if(j==length(temp.map) || all(is.na(temp.map[-(1:j)]))) {
          z <- max((1:length(temp.map))[-o])
          temp.map[j] <- max(temp.map,na.rm=TRUE)+(j-z+1)
        }
        else {
          temp.map[j] <- (min(temp.map[-(1:j)],na.rm=TRUE)+
                          max(temp.map[1:(j-1)],na.rm=TRUE))/2
        }
      }
    }
    
    names(temp.map) <- mnames[chr==uchr[i]]

    # pull out appropriate portion of genotype data
    data <- allgeno[,min.mar:(length(temp.map)+min.mar-1)]
    min.mar <- min.mar + length(temp.map)
    colnames(data) <- names(temp.map)

    geno[[i]] <- list(data=data,map=temp.map)
    if(uchr[i] == "X" || uchr[i] == "x") 
      class(geno[[i]]) <- "X"
    else class(geno[[i]]) <- "A"
  }
  colnames(pheno) <- pnames

  # check that data dimensions match 
  n.mar1 <- sapply(geno,function(a) ncol(a$data))
  n.mar2 <- sapply(geno,function(a) length(a$map))
  n.phe <- ncol(pheno)
  n.ind1 <- nrow(pheno)
  n.ind2 <- sapply(geno,function(a) nrow(a$data))
  if(any(n.ind1 != n.ind2)) {
    print(c(n.ind1,n.ind2))
    stop("Number of individuals in genotypes and phenotypes do not match.");
  }
  if(any(n.mar1 != n.mar2)) {
    print(c(n.mar,n.mar2))
    stop("Numbers of markers in genotypes and marker names files do not match.");
  }

  # print some information about the amount of data read
  cat(" --Read the following data:\n");
  cat("\t", n.ind1, " individuals\n");
  cat("\t", sum(n.mar1), " markers\n");
  cat("\t", n.phe, " phenotypes\n");

  # determine map type: f2 or bc or 4way?
  if(max(allgeno[!is.na(allgeno)])<=2) type <- "bc"
  else type <- "f2"
  cross <- list(geno=geno,pheno=pheno)
  class(cross) <- c(type,"cross")

  # check that nothing is strange in the genotype data
  cross.type <- class(cross)[1]
  if(cross.type=="f2") max.gen <- 5
  else max.gen <- 2

  u <- unique(allgeno)
  if(any(!is.na(u) & (u > max.gen | u < 1))) 
    stop(paste("There are stange values in the genotype data :",
               paste(sort(u),collapse=":"), "."))

  cross
}
  

######################################################################
#
# read.cross.csv
#
# read data in comma-delimited format
#
######################################################################

read.cross.csv <-
function(dir,file,sep=",",na.strings="-",genotypes=c("A","H","B","C","D"))
{
  # create file names 
  if(missing(file)) file <- "data.csv"

  if(!missing(dir)) {
    # remove ending "/" if it exists
    n <- nchar(dir)
    if(substr(dir,n,n) == "/")
      dir <- substr(dir,0,n-1)
    file <- paste(dir,file, sep="/")
  }

  # read data
  data <- as.matrix(read.table(file,sep=sep,na.strings=na.strings,
                               as.is=TRUE))
  data[!is.na(data) & data == ""] <- NA

  # pull apart phenotypes, genotypes and map
  n.phe <- max((1:ncol(data))[is.na(data[2,])])
  pheno <- matrix(as.numeric(data[-(1:3),1:n.phe]),nrow=nrow(data)-3)
  colnames(pheno) <- data[1,1:n.phe]
  mnames <- data[1,-(1:n.phe)]
  chr <- data[2,-(1:n.phe)]
  map <- as.numeric(data[3,-(1:n.phe)])
  allgeno <- matrix(match(data[-(1:3),-(1:n.phe)],genotypes),
                    ncol=ncol(data)-n.phe,nrow=nrow(data)-3)

  # fix up map information
  # number of chromosomes
  uchr <- unique(chr)
  n.chr <- length(uchr)
  geno <- vector("list",n.chr)
  names(geno) <- uchr
  min.mar <- 1
  for(i in 1:n.chr) { # loop over chromosomes
    # create map
    temp.map <- map[chr==uchr[i]]
    names(temp.map) <- mnames[chr==uchr[i]]

    # pull out appropriate portion of genotype data
    data <- allgeno[,min.mar:(length(temp.map)+min.mar-1)]
    min.mar <- min.mar + length(temp.map)
    colnames(data) <- names(temp.map)

    geno[[i]] <- list(data=data,map=temp.map)
    if(uchr[i] == "X" || uchr[i] == "x") 
      class(geno[[i]]) <- "X"
    else class(geno[[i]]) <- "A"
  }

  # check that data dimensions match 
  n.mar1 <- sapply(geno,function(a) ncol(a$data))
  n.mar2 <- sapply(geno,function(a) length(a$map))
  n.phe <- ncol(pheno)
  n.ind1 <- nrow(pheno)
  n.ind2 <- sapply(geno,function(a) nrow(a$data))
  if(any(n.ind1 != n.ind2)) {
    print(c(n.ind1,n.ind2))
    stop("Number of individuals in genotypes and phenotypes do not match.");
  }
  if(any(n.mar1 != n.mar2)) {
    print(c(n.mar,n.mar2))
    stop("Numbers of markers in genotypes and marker names files do not match.");
  }

  # print some information about the amount of data read
  cat(" --Read the following data:\n");
  cat("\t", n.ind1, " individuals\n");
  cat("\t", sum(n.mar1), " markers\n");
  cat("\t", n.phe, " phenotypes\n");

  # determine map type: f2 or bc or 4way?
  if(max(allgeno[!is.na(allgeno)])<=2) type <- "bc"
  else type <- "f2"
  cross <- list(geno=geno,pheno=pheno)
  class(cross) <- c(type,"cross")

  # check that nothing is strange in the genotype data
  cross.type <- class(cross)[1]
  if(cross.type=="f2") max.gen <- 5
  else max.gen <- 2

  u <- unique(allgeno)
  if(any(!is.na(u) & (u > max.gen | u < 1))) 
    stop(paste("There are stange values in the genotype data :",
               paste(sort(u),collapse=":"), "."))

  cross
}
  

# end of read.cross.R
