#####################################################################
#
# qtlcart_io.R
#
# copyright (c) 2002-3, Brian S. Yandell
#          [modified by Karl W. Broman and Hao Wu]
# last modified Jun, 2003
# first written Jun, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/qtl package
# Contains: read.cross.qtlcart, read.cro.qtlcart, read.map.qtlcart,
#           write.cross.qtlcart
##############################################################################

read.cross.qtlcart <-
function (dir, crofile, mapfile)
{
    if (missing(mapfile)) stop("Missing mapfile.")
    if (missing(crofile)) stop("Missing crofile.")

    if(!missing(dir) && dir != "") {
      mapfile <- file.path(dir, mapfile)
      crofile <- file.path(dir, crofile)
    }
    map <- read.map.qtlcart( mapfile )
    cro <- read.cro.qtlcart( crofile )

    cat(" --Read the following data:\n")
    cat("       Type of cross:         ", cro$cross, "\n")
    cat("       Number of individuals: ", nrow( cro$markers ), "\n")
    cat("       Number of markers:     ", ncol( cro$markers ), "\n")
    cat("       Number of phenotypes:  ", nrow( cro$traits ), "\n")

    maplen <- unlist(lapply(map,length))
    markers <- split( as.data.frame( t( cro$markers )),
                     ordered( rep(names( maplen ), maplen )))

    Geno <- list()
    for( i in names( map )) {
      name.markers <- names( map[[i]] )
      markers[[i]] <- t( markers[[i]] )
      colnames( markers[[i]] ) <- name.markers
      tmp <- list( data = markers[[i]], map = map[[i]] )
#      class( tmp ) <- if( length( grep( "^.*[xX].*", name.markers )))
#        "X"
#      else
#        "A"
      # determine whether autosomal chromosome or X chromosome
      #     using the chromosome name
      class(tmp) <- ifelse(length(grep("[Xx]", i)), "X", "A")
      Geno[[i]] <- tmp
    }
    cross <- list(geno = Geno, pheno = cro$traits )
    class(cross) <- c( cro$cross, "cross")

    cross$pheno <- as.data.frame(cross$pheno)

    list(cross,FALSE)
}

read.cro.qtlcart <-
function( file )
{
  ## translation from cro to R/qtl (see read.cross)
  ## -1	NA	missing data
  ##  0	1	AA
  ##  1	2	AB
  ##  2	3	BB
  ## 10	4	AA or AB
  ## 12	5	AB or BB
  ##
  f <- scan( file, what = "", blank.lines.skip = FALSE, sep = "\n", quiet = TRUE )
  ctrl <- seq( f )[ substring( f, 1, 1 ) == "-" ]
  s <- strsplit( f[ctrl], " " )
  ns <- character( length( ctrl ))
  for( i in seq( ctrl )) {
    ns[i] <- substring( s[[i]][1], 2 )
    s[[i]] <- s[[i]][ "" != s[[i]] ][-1]
  }
  names( s ) <- ns
  size <- as.numeric( s$n[1] )
  nmarkers <- as.numeric( s$p[1] ) - 1
  ntraits <- as.numeric( s$traits[1] )

  # cross type
  fix.ridh <- FALSE # indicator of whether to fix genotypes
  cross <- s$cross[1]
  if(cross=="RI1") {
    cross <- "riself"
    fix.ridh <- TRUE
  }
  else if(cross=="RI2") {
    cross <- "risib"
    fix.ridh <- TRUE
  }
  else if(cross=="RI0") {
    cross <- "bc" # doubled haploid
    fix.ridh <- TRUE
  }
  else if(cross=="B1" || cross=="B2") cross <- "bc"
  else if(cross=="SF2" || cross=="RF2") cross <- "f2"
  else if(cross!="f2" && cross!="bc" && cross!="f2ss" &&
          cross!="risib" && cross!="riself" && cross!="4way") {
    err <- paste("Cross type",cross,"not supported.")
    stop(err)
  }

  notraits <- as.numeric( s$otraits[1] )
  skip <- ctrl[ "s" == ns ]
  nlines <- ctrl[ "e" == ns ] - skip - 1
  trait.names <- f[ ctrl[ "Names" == ns ] + 1:ntraits ]
  ns <- strsplit( trait.names, " " )
  for( i in seq( ns ))
    ns[[i]] <- ns[[i]][ length( ns[[i]] ) ]
  trait.names <- unlist( ns )
  f <- matrix( scan( file, skip = skip, nlines = nlines, na.strings = ".",
                    blank.lines.skip = FALSE, quiet = TRUE ),
              ncol = size )

  traits <- t( f[-(1:(2+nmarkers)),] )
  if( nrow( traits ) == 1 )
    traits <- t( traits )
  dimnames( traits ) <- list( NULL, trait.names )
  f <- t( f[ 3:(2+nmarkers), ] )
  ## here is the translation
  f[ !is.na( f ) ] <- c(NA,1:3,rep(NA,7),4,NA,5)[ 2 + f[ !is.na( f ) ] ]

  if(fix.ridh && all(is.na(f) || f==1 || f==3))
    f[!is.na(f) & f==3] <- 2

  list( traits = traits, markers = f, cross = cross )
}

read.map.qtlcart <-
function( file )
{
## only interested in chromosomes, marker IDs and positions
  f <- scan( file, what = "", blank.lines.skip = FALSE, sep = "\n", quiet = TRUE )
  ctrl <- seq( f )[ substring( f, 1, 1 ) == "-" ]
  getvalue <- function( s, f, ctrl ) {
    tmp <- unlist( strsplit( f[ ctrl[ substring( f[ctrl], 2, 3 ) == s ] ],
                            " " ))
    as.numeric( tmp[ "" != tmp ][2] )
  }
  nchrom <- getvalue( "c ", f, ctrl )
  nmarkers <- getvalue( "i ", f, ctrl )

  ## marker positions
  tmp <- range( seq( f )[ substring( f, 1, 3 ) == "-l " ] )
  s <- strsplit( f[ tmp[1] ], "" )[[1]]
  b <- grep( "|", s, extended = FALSE )
  s <- grep( "0", s )
  s <- ceiling(( s[ length( s ) ] - s[2] ) / nchrom )
  position <- as.matrix( read.fwf( file, c( 1 + b, rep( s, nchrom )),
                                  skip = tmp[1]-1, n = tmp[2] )[,-1] )
  tmp <- grep( "-b", f )
  markers <- scan( file, list(1,2,""), skip = tmp[1], nlines = nmarkers,
                  blank.lines.skip = FALSE, quiet = TRUE )
  chroms <- scan( file, list(1,""), skip = tmp[2], nlines = nchrom,
                 blank.lines.skip = FALSE, quiet = TRUE )[[2]]

  map <- list( )
  for( i in seq( nchrom )) {
    tmp <- cumsum( position[ !is.na( position[,i] ), i ] )
    tmp <- tmp[ - length( tmp ) ]
    names( tmp ) <- markers[[3]][ i == markers[[1]] ]
    map[[ chroms[i] ]] <- tmp
  }
  map
}

write.cross.qtlcart <-
function( cross, filestem="data", chr )
{
  require( qtl )
  if(!missing(chr))
    cross <- subset(cross,chr=chr)

  n.ind <- nind(cross)
  tot.mar <- totmar(cross)
  n.phe <- nphe(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)

  type <- class(cross)[1]
  if(type=="bc") type <- "B1"
  else if(type=="f2" || type=="f2ss") type <- "RF2"
  else if(type=="riself") type <- "RI1"
  else if(type=="risib") type <- "RI2"
  else {
    warn <- paste("Cross type", type, "may not work with QTL Cartographer.")
    warning(warn)
  }

  # write genotype and phenotype data
  file <- paste(filestem, ".cro", sep="")
  if( file.exists( file )) {
    warning( paste( "previous file", file, "moved to *.mov" ))
    file.rename( file, paste( file, "mov", sep = "." ))
  }
  write("#  123456789 -filetype Rcross.out", file, append=FALSE)

  ## write numbers of progeny, markers and phenotypes
  write( paste( "-n   ", n.ind ), file, append=TRUE)
  write( paste( "-p   ", 1 + tot.mar ), file, append=TRUE)
  ## write experiment type
  write( paste( "-cross", type ), file, append=TRUE)

  ## write numbers of progeny, markers and phenotypes
  write( paste( "-traits   ", n.phe ), file, append=TRUE)
  write( "-Names of traits...", file, append=TRUE)
  phe <- names( cross$pheno )
  for( i in seq( phe ))
    write( paste( i, phe[i] ), file, append=TRUE)
  write( paste( "-otraits   ", 0 ), file, append=TRUE)

  ## write genotype and phenotype data by individual
  write( "-s", file, append=TRUE)
  for( ind in 1:n.ind ) {
    write( paste( ind, 1 ), file, append=TRUE)
    for(i in 1:n.chr) {
      g <- unlist( cross$geno[[i]]$data[ind,] )
      g[ is.na( g ) ] <- 0
      g <- c(-1,0,1,2,10,12)[ 1 + g ]
#      if( length( g ) < 40) # replaced by Hao with the following
      if( length( g ) <= 40)
        write(paste( "      ", paste( g, collapse = " " )), file, append=TRUE)
      else {
#        lo <- seq( 1, n.ind-1, by=40) # replaced by Hao with the following
        lo <- seq( 1, length(g), by=40)
#        hi <- c( lo[-1]+1, length( g )) # replaced by Hao with the following
        hi <- c( lo[-1]-1, length( g ))
        for(k in seq(along=lo)) {
          write( paste( "      ", paste( g[lo[k]:hi[k]], collapse = " " )),
                file, append=TRUE)
        }
      }
    } ## end writing marker data
    p <- c( cross$pheno[ind,])
    tmp <- format( p )
    tmp[ is.na( p ) ] <- "."
    write( paste( "       ", tmp ), file, append = TRUE )
    ## end of writing phenotype data
  }
  write( "-e", file, append = TRUE )
  write( "-q", file, append = TRUE )
#  unlink( file )

  # make "prep" file with map information
  file <- paste(filestem, ".map", sep="")
  if( file.exists( file )) {
    warning( paste( "previous file", file, "moved to *.mov" ))
    file.rename( file, paste( file, "mov", sep = "." ))
  }
  write("#  123456789 -filetype Rmap.out", file, append=FALSE)

  ## write numbers of progeny, markers and phenotypes
  write( "-s", file, append=TRUE)
  write( "-f 1", file, append=TRUE)
  write( "-p 0.0000", file, append=TRUE)
  write( "-u c", file, append=TRUE)
  write( "#", file, append=TRUE)

  write( paste( "-c", n.chr ), file, append=TRUE)
  write( paste( "-i", tot.mar ), file, append=TRUE)

  map <- lapply( cross$geno, function( x ) x$map )
  maplen <- unlist( lapply( map, length ))
  ## mean and SD of number of markers
  write( paste( "-m", round( mean( maplen ), 3 )), file, append=TRUE)
  write( paste( "-vm", round( sqrt( var( maplen )), 3 )), file, append=TRUE)

  mapdif <- lapply( map, diff )
  ## mean and SD of intermarker distances
  write( paste( "-d", round( mean( unlist( mapdif )), 3 )), file, append=TRUE)
  write( paste( "-vd", round( sqrt( var( unlist( mapdif ))), 3 )), file, append=TRUE)
  write( "-t 0.0000", file, append=TRUE)
  write( "#", file, append=TRUE)
  write( "          |   Chromosome----> ", file, append=TRUE)
  write( "--------------------------------------", file, append=TRUE)
  mapmat <- matrix( NA, 1 + max( maplen ), n.chr )
  mapmat[ 1, ] <- 0
  for( i in seq( along = maplen )) {
    tmp <- c( mapdif[[i]],0)
    mapmat[1 + seq( along = tmp ), i ] <- tmp
  }
  mapmat <- format( mapmat )
  ncmap <- nchar( mapmat[1] )
  mapmat[ grep( "NA", mapmat ) ] <- paste( rep( " ", ncmap ), collapse = "" )
  tmp <- format( seq( n.chr ))
  write( paste( "Marker    | ",
               paste( tmp, collapse =
                     paste( rep( " ", max( 1, ncmap - nchar( tmp ))), collapse = "" ))),
        file, append=TRUE)
  write( "--------------------------------------", file, append=TRUE)
  for( i in seq( nrow( mapmat )))
    write( paste( "-l     ", i - 1, "|",
                 paste( mapmat[i,], collapse = " " )),
          file, append=TRUE)

  write( "---------------------------------------", file, append=TRUE)
  write( paste( "-Number   |", paste( maplen, collapse = "  " )),
        file, append=TRUE)

  write( "Names and positions of the markers", file, append=TRUE)
  write( "Chrom  Mark  Name", file, append=TRUE)
  write( "-b MarkerNames", file, append=TRUE)
  for( i in 1:n.chr )
    for( j in seq( along = map[[i]] ))
      write( paste( i, j, names( map[[i]] )[j] ), file, append=TRUE)
  write( "-e MarkerNames", file, append=TRUE)
  write( "Names of the Chromosomes", file, append=TRUE)
  write( "-b ChromosomeNames", file, append=TRUE)
  for( i in 1:n.chr )
    write( paste( i, names( map )[i] ), file, append=TRUE)
  write( "-e ChromosomeNames", file, append=TRUE)
}

# end of qtlcart_io.R
