######################################################################
#
# summary.cross.R
#
# copyright (c) 2001, Karl W Broman, Johns Hopkins University
# last modified Nov, 2001
# first written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: summary.cross, print.summary.cross, nind, nchr, nmar,
#           totmar, nphe
#
######################################################################

summary.cross <-
function(object,...)
{
#  if(is.na(match("cross",class(object))))
#    stop("This is not an object of class cross.")
    
  n.ind <- nind(object)
  tot.mar <- totmar(object)
  n.phe <- nphe(object)
  n.chr <- nchr(object)
  n.mar <- nmar(object)
  type <- class(object)[1]

  Geno <- object$geno[[1]]$data
  if(n.chr > 1)
    for(i in 2:n.chr)
      Geno <- cbind(Geno,object$geno[[i]]$data)

  missing.gen <- mean(is.na(Geno))
  
  if(type=="f2") {
    typings <- table(factor(Geno[!is.na(Geno)], levels=1:5))
    names(typings) <- c("AA","AB","BB","not BB","not AA")
  }
  else if(type=="bc") {
    typings <- table(factor(Geno[!is.na(Geno)], levels=1:2))
    names(typings) <- c("AA","AB")
  }
  else 
    typings <- table(factor(Geno[!is.na(Geno)]))

  typings <- typings/sum(typings)

  missing.phe <- as.numeric(cbind(apply(object$pheno,2,function(a) mean(is.na(a)))))

  # check that object$geno[[i]]$data has colnames and that they match
  #     the names in object$geno[[i]]$map
  for(i in 1:n.chr) {
    nam1 <- colnames(object$geno[[i]]$data)
    map <- object$geno[[i]]$map
    if(is.matrix(map)) nam2 <- colnames(map)
    else nam2 <- names(map)
    chr <- names(object$geno)[[i]]
    if(is.null(nam1)) {
      warn <- paste("The data matrix for chr", chr,
                    "lacks column names")
      warning(warn)
    }
    if(is.null(nam2)) {
      warn <- paste("The genetic map for chr", chr,
                    "lacks column names")
      warning(warn)
    }
    if(any(nam1 != nam2)) {
      warn <- paste("Marker names in the data matrix and genetic map\n",
                    "for chr ", chr, " do not match.",sep="")
      stop(warn)
    }
      
    if((is.matrix(map) && (any(diff(map[1,])<0) || any(diff(map[2,])<0))) ||
       (!is.matrix(map) && any(diff(map)<0)))
        stop(paste("Markers out of order on chr", chr))
  }
    
  if(!is.data.frame(object$pheno))
    warning("Phenotypes should be a data.frame.")

  cross.summary <- list(type=type, n.ind = n.ind, n.phe=n.phe, 
			n.chr=n.chr, n.mar=n.mar,
			missing.gen=missing.gen,typing.freq=typings,
			missing.phe=missing.phe)
  class(cross.summary) <- "summary.cross"
  cross.summary
  
}


print.summary.cross <-
function(x,...)
{
  cat("\n")
  if(x$type=="f2") cat("    F2 intercross\n\n")
  else if(x$type=="bc") cat("    Backcross\n\n")
  else if(x$type=="4way") cat("    4-way cross\n\n")
  else cat(paste("    cross", x$type, "\n\n",sep=" "))

  cat("    No. individuals: ", x$n.ind,"\n\n")
  cat("    No. phenotypes:  ", x$n.phe,"\n\n")
  cat("    No. chromosomes: ", x$n.chr,"\n")
  cat("    Total markers:   ", sum(x$n.mar), "\n")
  cat("    No. markers:     ", x$n.mar, "\n")
  cat("\n")
  cat("    Percent genotyped: ", round((1-x$missing.gen)*100,1), "\n")
  cat("    Genotypes (%):     ", 
      paste(names(x$typing.freq),round(x$typing.freq*100,1),sep=":", collapse="  "),
      "\n")
  cat("\n")
  cat("    Percent phenotyped: ", round((1-x$missing.phe)*100,1), "\n")
  cat("\n")
}



nind <-
function(object)
{
  if(any(is.na(match(c("pheno","geno"),names(object)))))
    stop("This is not an object of class cross.")

  n.ind1 <- nrow(object$pheno)
  n.ind2 <- sapply(object$geno,function(x) nrow(x$data))
  if(any(n.ind2 != n.ind1))
    stop("Different numbers of individuals in genotypes and phenotypes.")
  n.ind1
}

nchr <-
function(object)
{
  if(any(is.na(match(c("pheno","geno"),names(object)))))
    stop("This is not an object of class cross.")
  length(object$geno)
}

nmar <- 
function(object)
{
  if(any(is.na(match(c("pheno","geno"),names(object)))))
    stop("This is not an object of class cross.")

  if(!is.matrix(object$geno[[1]]$map))
    n.mar1 <- sapply(object$geno, function(x) length(x$map))
  else # sex-specific maps
    n.mar1 <- sapply(object$geno, function(x) ncol(x$map))

  n.mar2 <- sapply(object$geno, function(x) ncol(x$data))
  if(any(n.mar1 != n.mar2))
    stop("Different numbers of markers in genotypes and maps.")
  n.mar1
}

totmar <-
function(object)
{
  if(any(is.na(match(c("pheno","geno"),names(object)))))
    stop("This is not an object of class cross.")

  if(!is.matrix(object$geno[[1]]$map))
    totmar1 <- sum(sapply(object$geno, function(x) length(x$map)))
  else # sex-specific maps
    totmar1 <- sum(sapply(object$geno, function(x) ncol(x$map)))
  totmar2 <- sum(sapply(object$geno, function(x) ncol(x$data)))
  if(totmar1 != totmar2)
    stop("Different numbers of markers in genotypes and maps.")
  totmar1
}

nphe <-
function(object) {
  if(any(is.na(match(c("pheno","geno"),names(object)))))
    stop("This is not an object of class cross.")

  ncol(object$pheno)
}

# end of summary.cross.R
