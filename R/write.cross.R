######################################################################
#
# write.cross.R
#
# copyright (c) 2000-2001, Karl W Broman, Johns Hopkins University
# last modified Nov, 2001
# first written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/qtl package
# Contains: write.cross, write.cross.mm, write.cross.csv
#
######################################################################


######################################################################
#
# write.cross: Wrapper for the other write.cross functions
#
######################################################################

write.cross <-
function(cross, format=c("csv","mm"), filestem="data", chr, digits=5)
{
  match.arg <- match.arg(format)
  if(missing(chr)) chr <- 1:nchr(cross)

  if(format=="csv") write.cross.csv(cross,filestem,chr,digits)
  else write.cross.mm(cross,filestem,chr,digits)
}


######################################################################
#
# write.cross.mm: Write data for an experimental cross in Mapmaker
#                 format
#
#           creates two files: "raw" file with geno & pheno data
#                              "prep" file with map information
#
######################################################################

write.cross.mm <-
function(cross, filestem="data", chr, digits=5)
{
  if(!missing(chr)) 
    cross <- subset(cross,chr=chr)

  n.ind <- nind(cross)
  tot.mar <- totmar(cross)
  n.phe <- nphe(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)
  
  type <- class(cross)[1]
  if(type != "f2" && type != "bc")
    stop("write.cross.mm only works for intercross and backcross data.")

  # write genotype and phenotype data
  file <- paste(filestem, ".raw", sep="")
  
  # write experiment type
  if(type == "f2") 
    write("data type f2 intercross", file, append=FALSE)
  else 
    write("data type f2 backcross", file, append=FALSE)

  # write numbers of progeny, markers and phenotypes
  write(paste(n.ind, tot.mar, n.phe), file, append=TRUE)

  # max length of marker name
  mlmn <- max(nchar(unlist(lapply(cross$geno,function(a) colnames(a$data)))))+1

  # write marker data
  for(i in 1:n.chr) {
    for(j in 1:n.mar[i]) {
      mn <- paste("*", colnames(cross$geno[[i]]$data)[j], sep="")
      if(nchar(mn) < mlmn)
        mn <- paste(mn,paste(rep(" ", mlmn-nchar(mn)),collapse=""),sep="")
      g <- cross$geno[[i]]$data[,j]
      
      x <- rep("", n.ind)
      x[is.na(g)] <- "-"
      x[!is.na(g) & g==1] <- "A"
      x[!is.na(g) & g==2] <- "H"
      if(type == "f2") {
        x[!is.na(g) & g==3] <- "B"
        x[!is.na(g) & g==4] <- "D"
        x[!is.na(g) & g==5] <- "C"
      }

      if(n.ind < 60)
        write(paste(mn, paste(x,collapse="")), file, append=TRUE)
      else {
        lo <- seq(1,n.ind-1,by=60)
        hi <- c(lo[-1]-1,n.ind)
        for(k in seq(along=lo)) {
          if(k==1) write(paste(mn,paste(x[lo[k]:hi[k]],collapse="")),file,append=TRUE)
          else write(paste(paste(rep(" ", mlmn),collapse=""),
                           paste(x[lo[k]:hi[k]],collapse="")),file,append=TRUE)
        }
      }

    }
  } # end writing marker data

  # max length of phenotype name
  mlpn <- max(nchar(colnames(cross$pheno)))+1

  # write phenotypes
  for(i in 1:n.phe) {
    pn <- paste("*",colnames(cross$pheno)[i],sep="")
    if(nchar(pn) < mlpn)
      pn <- paste(pn, paste(rep(" ", mlpn-nchar(pn)),collapse=""),sep="")

    x <- as.character(round(cross$pheno[,i],digits))
    x[x=="NA"] <- "-"

    if(n.ind < 10)
      write(paste(pn, paste(x,collapse="")), file, append=TRUE)
    else {
      lo <- seq(1,n.ind-1,by=10)
      hi <- c(lo[-1]+1,n.ind)
      for(k in seq(along=lo)) {
        if(k==1) write(paste(pn,paste(x[lo[k]:hi[k]],collapse=" ")),file,append=TRUE)
        else write(paste(paste(rep(" ", mlpn),collapse=""),
                         paste(x[lo[k]:hi[k]],collapse=" ")),file,append=TRUE)
      }
    }
  }
    

  # make "prep" file with map information
  file <- paste(filestem, ".prep", sep="")

  for(i in 1:n.chr) {
    cname <- paste("chr", names(cross$geno)[i], sep="")
    line <- paste("make chromosome", cname)
    if(i==1) write(line, file, append=FALSE)
    else write(line, file, append=TRUE)

    mn <- names(cross$geno[[i]]$map)
    dis <- round(diff(cross$geno[[i]]$map),2)
    dis <- paste("=", dis, sep="")
    write(paste(paste("sequence", mn[1]), paste(dis,mn[-1],collapse=" ")),
          file, append=TRUE)

    write(paste("anchor", cname), file, append=TRUE)
    

  }

} 

######################################################################
#
# write.cross.csv: Write data for an experimental cross in
#                  comma-delimited format (the same format as is read
#                  by read.cross.csv)
#
######################################################################

write.cross.csv <-
function(cross, filestem="data", chr, digits=5)
{
  if(!missing(chr)) cross <- subset(cross,chr=chr)

  n.ind <- nind(cross)
  tot.mar <- totmar(cross)
  n.phe <- nphe(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)
  
  type <- class(cross)[1]
  if(type != "f2" && type != "bc")
    stop("write.cross.csv only works for intercross and backcross data.")

  file <- paste(filestem, ".csv", sep="")
  
  geno <- matrix(ncol=tot.mar,nrow=n.ind)
  alleles <- c("A","H","B","D","C")
  firstmar <- 1
  for(i in 1:n.chr) {
    # replace allele numbers with 
    geno[,firstmar:(firstmar+n.mar[i]-1)] <-
      alleles[match(cross$geno[[i]]$data,1:5)]
    firstmar <- firstmar + n.mar[i]
  }
  geno[geno=="NA"] <- "-"
  data <- cbind(matrix(as.character(round(unlist(cross$pheno),digits)),nrow=n.ind),geno)
  colnames(data) <- c(colnames(cross$pheno),
                      unlist(lapply(cross$geno, function(a) colnames(a$data))))
  chr <- rep(names(cross$geno),n.mar)
  pos <- unlist(lapply(cross$geno,function(a) a$map))
  chr <- c(rep("",n.phe),chr)
  pos <- c(rep("",n.phe),as.character(round(pos,digits)))

  # write names
  write.table(matrix(colnames(data),nrow=1),file,append=FALSE,quote=FALSE,
              sep=",",row.names=FALSE,col.names=FALSE)
  # write chr IDs
  write.table(matrix(chr,nrow=1),file,append=TRUE,quote=FALSE,sep=",",
              row.names=FALSE,col.names=FALSE)
  # write marker positions
  write.table(matrix(pos,nrow=1),file,append=TRUE,quote=FALSE,sep=",",
              row.names=FALSE,col.names=FALSE)
  # write phenotype and genotype data
  write.table(data,file,append=TRUE,quote=FALSE,sep=",",row.names=FALSE,
              col.names=FALSE)

}

# end of write.cross.R 
