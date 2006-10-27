######################################################################
#
# makeqtl.R
#
# copyright (c) 2002-6, Hao Wu, The Jackson Laboratory
#                     and Karl W. Broman, Johns Hopkins University
# last modified Oct, 2006
# first written Apr, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: makeqtl, replaceqtl, addqtl, dropqtl, locatemarker
#
######################################################################

######################################################################
#
# This is the function to construct an object of class "qtl"
# The phenotype data and genotype data for a given list of
# chromosome and locations will be extracted from the input
# "cross" object
#
######################################################################

######################################################################
#
# Notes/Question:
#  1. Do we want to pull out draws and prob at the same time?
#     If user specifed pos is far from the real marker, there
#     might be some problem for pulling out genoprob
#  2. The utility functions can be put in util.R
#
######################################################################

makeqtl <-
  function(cross, chr, pos, qtl.name)
{
  if( !sum(class(cross) == "cross") )
    stop("The first input variable must be  an object of class cross")

  # cross type
  type <- class(cross)[1]
  chrtype <- sapply(cross$geno, class)
  sexpgm <- getsex(cross)
  
  # check phenotypes
#  if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
#  if(pheno.col < 1 || pheno.col > nphe(cross))
#    stop("Specified phenotype column is invalid.")
  
  # chr, pos and qtl.name must have the same length
  if(length(chr) != length(pos))
    stop("Input chr and pos must have the same length.")
  else if( !missing(qtl.name) )
    if( length(chr) != length(qtl.name) )
      stop("Input chr and qtl.name must have the same length.")

  # local variables
  n.ind <- nrow(cross$pheno) # number of individuals
  n.pos <- length(chr) # number of selected markers
  n.gen <- NULL

  # initialize output object
  qtl <- NULL
  
  # take out the imputed genotypes and/or genoprobs for the
  # selected markers (if there are there)
  if("draws" %in% names(cross$geno[[1]])) { # draws is there
    # take out imputed genotype data
    n.draws <- dim(cross$geno[[1]]$draws)[3] # number of draws
    # initialize geno matrix for selected markers
    geno <- array(rep(0, n.ind*n.pos*n.draws),
                  dim=c(n.ind, n.pos, n.draws))
    for(i in 1:n.pos) {
      # get the index for this chromosome
      i.chr <- which(chr[i]==names(cross$geno))
      if(length(i.chr) == 0) # no this chromosome in cross 
        stop("There's no chromosome number ", chr[i], " in input cross object")
      i.pos <- pos[i] # marker position

      # make the genetic map for this chromosome
      if("map" %in% names(attributes(cross$geno[[i.chr]]$draws)))
        map <- attr(cross$geno[[i.chr]]$draws,"map")
      else {
        stp <- attr(cross$geno[[i.chr]]$draws, "step")
        oe <- attr(cross$geno[[i.chr]]$draws, "off.end")
      
        if("stepwidth" %in% names(attributes(cross$geno[[i.chr]]$draws)))
          stpw <- attr(cross$geno[[i.chr]]$draws, "stepwidth")
        else stpw <- "fixed"
        map <- create.map(cross$geno[[i.chr]]$map,stp,oe,stpw)
      }

      # pull out the female map if there are sex-specific maps
      if(is.matrix(map)) map <- map[1,]

      # locate this marker (given chromosome and position)
      marker.idx <- locatemarker(map, i.pos, i.chr, flag="draws")
      
      # if everything is all right, take the genotype
      geno[,i,] <- cross$geno[[i.chr]]$draws[,marker.idx,]

      # no. genotypes
      n.gen[i] <- length(getgenonames(type,chrtype[i.chr],"full",sexpgm, attributes(cross)))
      
      # Fix up X chromsome here
      if(chrtype[i.chr]=="X")
        geno[,i,] <- reviseXdata(type,"full",sexpgm,draws=geno[,i,,drop=FALSE],
                                 cross.attr=attributes(cross))
    }
    # give geno dimension names
    # the 2nd dimension called "Q1", "Q2", etc.
    dimnames(geno) <- list(NULL, paste("Q", 1:n.pos, sep=""), NULL)
    # output 
    qtl$geno <- geno
  }
  
  if("prob" %in% names(cross$geno[[1]])) { # prob is there
    # initialize prob matrix
    prob <- vector("list",n.pos)

    # locate the marker
    for(i in 1:n.pos) {

      # get the index for this chromosome
      i.chr <- which(chr[i]==names(cross$geno))
      if(length(i.chr) == 0) # no this chromosome in cross
        stop("There's no chromosome number ", chr[i], " in input cross object")
      i.pos <- pos[i] # marker position

      # locate this marker (given chromosome and position)
      marker.idx <- locatemarker(cross$geno[[i.chr]]$map, i.pos, i.chr, flag="prob")

      # take genoprob
      prob[[i]] <- cross$geno[[i.chr]]$prob[,marker.idx,]

      ### Fix up X chromsome here ###
    }
    qtl$prob <- prob
  }

  if( sum(c("draws","prob") %in% names(cross$geno[[1]]))==0 ) 
    stop("You need to run calc.genoprob() or sim.geno() first.")

  if(missing(qtl.name))  # no given qtl names
    # make qtl names
    qtl.name <- paste( paste("Chr",chr,sep=""), pos, sep="@")

  # output object
  qtl$name <- qtl.name
  qtl$chr <- chr
  qtl$pos <- pos
  qtl$n.qtl <- n.pos
  qtl$n.ind <- nind(cross)
  qtl$n.gen <- n.gen

  class(qtl) <- "qtl"
  
  qtl
}
  


######################################################################
#
# This is the function to replace one QTL by another.
# This is the internal function and not supposed to be used by user
#
######################################################################
replaceqtl <-
  function(cross, qtl, replace, by.chr, by.pos, by.name, map)
{
  # update QTL name
  if(missing(by.name))
    by.name <- paste(paste("Chr",by.chr,sep=""), by.pos, sep="@")
  qtl$name[replace] <- by.name

  # update chr and pos
  qtl$chr[replace] <- by.chr
  qtl$pos[replace] <- by.pos
  
  # update the imputed genotype and n.gen vector (if any)
  if("geno" %in% names(qtl)) {
    if(missing(map))  { # make genetic map on this chromosome

      if("map" %in% names(attributes(cross$geno[[by.chr]]$draws)))
        map <- attr(cross$geno[[by.chr]]$draws,"map")
      else {
        stp <- attr(cross$geno[[by.chr]]$draws, "step")
        oe <- attr(cross$geno[[by.chr]]$draws, "off.end")
      
        if("stepwidth" %in% names(attributes(cross$geno[[by.chr]]$draws)))
          stpw <- attr(cross$geno[[by.chr]]$draws, "stepwidth")
        else stpw <- "fixed"
        map <- create.map(cross$geno[[by.chr]]$map,stp,oe,stpw)
      }

      # pull out female map in case that there are sex-specific maps
      if(is.matrix(map)) map <- map[1,]
    }

    # locate this marker (given chromosome and position)
    marker.idx <- locatemarker(map, by.pos, by.chr, "draws")

    # replace the genotypes
    qtl$geno[,replace,] <- cross$geno[[by.chr]]$draws[,marker.idx,]

     # update number of genotypes
    type <- class(cross)[[1]]
    if(type == "f2") {
      if(class(cross$geno[[by.chr]]) == "A") # autosomal
        qtl$n.gen[replace] <- 3
      else                             # X chromsome 
        qtl$n.gen[replace] <- 2
    }
    else if(type == "bc" || type=="risib" || type=="riself") 
      qtl$n.gen[replace] <- 2
    else if(type == "4way") 
      qtl$n.gen[replace] <- 4
    else 
      stop("replaceqtl not available for cross ", type)
  }
  
  # update the genoprob (if any)
  if("prob" %in% names(qtl)) {
    #locate the marker
    marker.idx <- locatemarker(cross$geno[[by.chr]]$map,
                                by.pos, by.chr, "prob")
    # replace genoprob
    qtl$prob[[replace]] <- cross$geno[[by.chr]]$prob[,marker.idx,]
  }

  # done
  
  qtl
}


######################################################################
#
# This is the function to add a QTL to given qtl object
# This is the internal function and not supposed to be used by user
#
######################################################################

addqtl <-
  function(cross, qtl, add.chr, add.pos, add.name, map)
{
  # update number of QTLs
  qtl$n.qtl <- qtl$n.qtl + 1

  # update chr and pos
  qtl$chr <- c(qtl$chr, add.chr)
  qtl$pos <- c(qtl$pos, add.pos)
  
  # add QTL name
  if(missing(add.name))
    add.name <- paste(paste("Chr",add.chr,sep=""), add.pos, sep="@")
  qtl$name[qtl$n.qtl] <- add.name
  
  # add new entry to the imputed genotype and n.gen vector (if any)
  if("geno" %in% names(qtl)) {
    # update number of genotypes
    type <- class(cross)[[1]]
    if(type == "f2") {
      if(class(cross$geno[[add.chr]]) == "A") # autosomal
        n.gen <- 3
      else                             # X chromsome 
        n.gen <- 2
    }
    else if(type == "bc" || type=="risib" || type=="riself") 
      n.gen <- 2
    else if(type == "4way") 
      n.gen <- 4
    else 
      stop("addqtl not available for cross ", type)
    qtl$n.gen <- c(qtl$n.gen, n.gen)
  
    # add the imputed genotype
    if(missing(map)) { # make genetic map on this chromosome, if missing

      if("map" %in% names(attributes(cross$geno[[add.chr]]$draws)))
        map <- attr(cross$geno[[add.chr]]$draws,"map")
      else {
        stp <- attr(cross$geno[[add.chr]]$draws, "step")
        oe <- attr(cross$geno[[add.chr]]$draws, "off.end")
      
        if("stepwidth" %in% names(attributes(cross$geno[[add.chr]]$draws)))
          stpw <- attr(cross$geno[[add.chr]]$draws, "stepwidth")
        else stpw <- "fixed"
        map <- create.map(cross$geno[[add.chr]]$map,stp,oe,stpw)
      }

      # pull out female map in case that there are sex-specific maps
      if(is.matrix(map)) map <- map[1,]
    }

    # locate this marker (given chromosome and position)
    marker.idx <- locatemarker(map, add.pos, add.chr, "draws")

    # reallocate memory for geno array
    n.ind <- dim(qtl$geno)[1]
    n.draw <- dim(qtl$geno)[3]
    geno <- array( rep(0, n.ind*n.draw*qtl$n.qtl),
                  c(n.ind, qtl$n.qtl, n.draw) )
  
    geno[,1:(qtl$n.qtl-1),] <- qtl$geno
    geno[,qtl$n.qtl,] <- cross$geno[[add.chr]]$draws[,marker.idx,]
    dimnames(geno) <- list(NULL, paste("Q", 1:qtl$n.qtl, sep=""), NULL)

    # replace geno in qtl
    qtl$geno <- geno
  }

  # add new entry to prob (if any)
  if("prob" %in% names(qtl)) {
    marker.idx <- locatemarker(cross$geno[[add.chr]]$map,
                                add.pos, add.chr, "prob")
    # reallocate memory for prob array
    n.ind <- dim(qtl$prob)[1]
    ngen <- dim(qtl$prob)[3]
    prob <- array( rep(0, n.ind*ngen*qtl$n.qtl),
                  c(n.ind, qtl$n.qtl, ngen))
    prob[,1:(qtl$n.qtl-1),] <- qtl$prob
    prob[,qtl$n.qtl,] <- cross$geno[[add.chr]]$prob[,marker.idx,]

    # replace prob in qtl
    qtl$prob <- prob
  }

  # done
  qtl
}

######################################################################
#
# This is the function to drop a QTL for a given qtl object
# This is the internal function and not supposed to be used by user
#
######################################################################
dropqtl <-
  function(qtl, drop)
{
  # input drop is an integer index
  # get the index for exclusing drop QTL
  idx <- setdiff(1:qtl$n.qtl, drop)
  
  # result object
  result <- NULL
  result$name <- qtl$name[idx]
  result$chr <- qtl$chr[idx]
  result$pos <- qtl$pos[idx]
  result$n.qtl <- qtl$n.qtl - 1
  result$n.ind <- qtl$n.ind
  result$n.gen <- qtl$n.gen[idx]
  result$geno <- qtl$geno[,idx,]
  result$prob <- qtl$prob[,idx,]
#  result$type <- type ## is this necessary? (and where should "type" come from?
  dimnames(result$geno) <- list(NULL, paste("Q", 1:result$n.qtl, sep=""),
                                NULL)

  class(result) <- "qtl"

  result
}


##################################################################
#
# locate the marker on a genetic map. Choose the nearest
# one if there's no marker or pseudomarker one the given
# location
#
# This is the internal function and not supposed to be used by user
#
###################################################################

locatemarker <-
  function(map, pos, chr, flag)
{  
  marker.idx <- which(map == pos)
  if( length(marker.idx)==0 ) {
    # there's no this marker, take the nearest marker instead
    # if there's a tie, take the first nearst one
    m.tmp <- abs(pos-map)
    marker.idx <- which(m.tmp==min(m.tmp))[[1]]
#    if(flag == "draws") {
#      msg <- "For draws: "
#    }
#    else if(flag == "prob") {
#      msg <- "For prob: "
#    }
#    msg <- paste(msg, "there's no marker on Chr ", chr, ", at ",
#                 pos,"cM.", sep="")
#    msg <- paste(msg, " Take marker at ", map[[marker.idx]], "cM instead.",
#                 sep="")
#    warning(msg)
  }

  if(length(marker.idx) > 1)
    marker.idx <- marker.idx[sample(length(marker.idx))]
  marker.idx
}

# end of makeqtl.R
