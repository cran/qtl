######################################################################
#
# fitqtl.R
#
# copyright (c) 2002, Hao Wu, The Jackson Laboratory
#                     and Karl W. Broman, Johns Hopkins University
# Last modified June, 2002
# first written April, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: fitqtl, parseformula, summary.fitqtl,
#           print.summary.fitqtl
#
######################################################################

######################################################################
#
# This is the function to fit a model and generate some tables
#
# Now only imputation method is implemented
#
#
######################################################################

fitqtl <-
function(pheno, qtl, covar=NULL, formula, method=c("imp"),
         dropone=TRUE)
{
  # some input checking stuff in here
  if( !sum(class(qtl) == "qtl") )
    stop("The first input variable must be  an object of class qtl")
  
  method <- match.arg(method)

  # check the input phenotypes and covarariates; drop individuals
  # with missing values.
  keep.ind <- !is.na(pheno)
  if(!is.null(covar)) { # loop thru covarariates
    for(i in 1:dim(covar)[2])
      keep.ind <- keep.ind & (!is.na(covar[,i]))
  }
  # if there IS missing data, do some subset
  if(any(!keep.ind)) {
    # subset pheno data
    pheno <- pheno[keep.ind]
    # subset covarariate
    covar.tmp <- as.data.frame(covar[keep.ind,])
    colnames(covar.tmp) <- colnames(covar)
    covar <- covar.tmp
    # hack input qtl object to drop individuals with missing data
    qtl$n.ind <- sum(keep.ind)
    qtl$geno <- qtl$geno[keep.ind,,]
  }
  
  # local variables
  n.ind <- qtl$n.ind # number of individuals
  n.qtl <- qtl$n.qtl # number of selected markers
  n.draws <- dim(qtl$geno)[3] # number of draws
  n.gen <- qtl$n.gen # number of genotypes
  
  if( is.null(covar) ){  # number of covarariates
    n.covar <- 0
  }
  else {
    n.covar <- dim(covar)[2]
  }
  
  # if formula is missing, build one
  # all QTLs and covarariates will be additive by default
  if(missing(formula)) {
    tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
    formula <- "y~Q1"
    if(n.qtl > 1) 
      for (i in 2:n.qtl) 
        formula <- paste(formula, tmp.Q[i], sep="+")
    if (n.covar) { # if covarariate is not empty
      tmp.C <- dimnames(covar)[[2]] # covarariate term names
      for(i in 1:n.covar)
        formula <- paste(formula, tmp.C[i], sep="+")
    }
    formula <- as.formula(formula)
  }

  # parse the input formula
  p <- parseformula(formula, dimnames(qtl$geno)[[2]], dimnames(covar)[[2]])

  # make an array n.gen.QC to represent the genotype numbers
  # for all input QTLs and covarariates. For covarariates the
  # number of genotyps is 1. This makes programming easier
  n.gen.QC <- c(n.gen[p$idx.qtl]-1, rep(1, p$n.covar))

  # covarariates to be passed to C function
  # This is done in case of that user input covar but has no covar in formula
  covar.C <- NULL
  if(!is.null(p$idx.covar))
    covar.C <- as.matrix(covar[,p$idx.covar])
  
  # call C function to do the genome scan
  if(method == "imp") {
    z <- .C("R_fitqtl_imp",
            as.integer(n.ind), # number of individuals
            as.integer(p$n.qtl), # number of qtls
            as.integer(n.gen.QC), # number of genotypes QTLs and covarariates
            as.integer(n.draws), # number of draws
            as.integer(qtl$geno[,p$idx.qtl,]), # genotypes for selected marker
            as.integer(p$n.covar), # number of covarariate
            as.double(covar.C), # covarariate
            as.integer(p$formula.intmtx),  # formula matrix for interactive terms
            as.integer(p$n.int), # number of interactions in the formula
            as.double(pheno), # phenotype
            # return variables
            lod=as.double(0), # LOD score
            df=as.integer(0), # degree of freedom

            PACKAGE="qtl")
  }

  ##### output ANOVA table for full model #####
  result.full <- matrix(NA, 3, 7)
  colnames(result.full) <- c("df", "SS", "MS", "LOD", "%var", "Pvalue(Chi2)",
                             "Pvalue(F)")
  rownames(result.full) <- c("Model", "Error", "Total")
  result.full[1,1] <- z$df # model degree of freedom
  # compute the SS for total
  Rss0 <- 0
  mpheno <- mean(pheno)
  for(i in 1:length(pheno)) {
    Rss0 <- Rss0 + (pheno[i]-mpheno)^2
  }
  # third row, for Total
  result.full[3,1] <- length(pheno) - 1 # total degree of freedom
  result.full[3,2] <- Rss0 # total sum of squares
    
  # first row, for Model
  result.full[1,1] <- z$df # df for Model
  # Variance explained by model
  result.full[1,5] <- 100 * (1 - exp(-2*z$lod*log(10)/n.ind))
  result.full[1,2] <- Rss0 * result.full[1,5]/100  # SS for model
  result.full[1,3] <- result.full[1,2]/z$df # MS for model
  result.full[1,4] <- z$lod # Model LOD score

  # Second row, for Error
  # df
  result.full[2,1] <- result.full[3,1] - result.full[1,1]
  # SS
  result.full[2,2] <- result.full[3,2] - result.full[1,2]
  # MS
  result.full[2,3] <- result.full[2,2] / result.full[2,1]

  # first row, P values
  # P value (chi2) for model
  result.full[1,6] <- 1 - pchisq(2*log(10)*z$lod, z$df)
  # P value (F statistics) for model
  df0 <- result.full[3,1]; df1 <- result.full[2,1];
  Rss1 <- result.full[2,2]
  Fstat <- ((Rss0-Rss1)/(df0-df1)) / (Rss1/df1)
  result.full[1,7] <- 1 - pf(Fstat, df0-df1, df1)

  ############# Finish ANOVA table for full model
  
  # initialize output object
  output <- NULL
  output$result.full <- result.full

  # drop one at a time?
  if(dropone) {
    # user wants to do drop one term at a time and output anova table

    # get the terms etc. for input formula
    f.terms <- terms(formula)
    f.order <- attr(f.terms, "order")
    f.label <- attr(f.terms, "term.labels")

    # initialize output matrix
    # ANOVA table will have five columns, e.g., df,Type III SS,
    # LOD, %var, Pvalue for each dropping term
    # Full model result will not be in this table
    result <- matrix(0, length(f.order), 6)
    colnames(result) <- c("df", "Type III SS", "LOD", "%var", "Pvalue(Chi2)",
                          "Pvalue(F)")
    rownames(result) <- rep("",length(f.order))

    # record the result for full model
#    result[1,1] <- z$df
#    result[1,3] <- z$lod
#    result[1,4] <-  100 * (1 - exp(-2*z$lod*log(10)/n.ind))
#    result[1,5] <- 1 - pchisq(2*log(10)*z$lod, z$df)
#    rownames(result)[1] <- "Full"
      
    drop.term.name <- NULL
    for( i in (1:length(f.order)) ) {
      # loop thru all terms in formula, from the highest orderp
      # the label of the term to be droped
      label.term.drop <- f.label[i]
      
      ### find the corresponding QTL name for this term ###
      # This is used for output ANOVA table
      if(f.order[i] == 1) {
        # this is a first order term
        # if the term label is like Q(q)1, Q(q)2, etc., then it's a QTL
        if( length(grep("Q[0-9]", label.term.drop, ignore.case=TRUE)) != 0) {
          idx.qtlname <- as.integer(substr(label.term.drop, 2, 10))
          drop.term.name[i] <- qtl$name[idx.qtlname]
        }
        else { # this is a covarariate
          drop.term.name[i] <- label.term.drop
        }
      }
      else {
        # this is a 2nd (or higher)order and the term is a string like "Q2:Q3:C1"
        # I use strsplit to split it to a string vector "Q2" "Q3" "C1".
        # then take out 2 and 3 as integer. Then find out the
        # QTL name from the input QTL object and concatenate them
        tmp.str <- strsplit(label.term.drop,":")[[1]]
        for(j in 1:length(tmp.str)) {
          if( length(grep("Q[0-9]", tmp.str[j], ignore.case=TRUE)) != 0 ) {
            # this is a QTL
            idx.qtlname <- as.integer(substr(tmp.str[j], 2, 100))
            tmp.str[j] <- qtl$name[idx.qtlname]
          }
          if(j == 1) # first term
            drop.term.name[i] <- tmp.str[j]
          else # not the first term
            drop.term.name[i] <- paste(drop.term.name[i], tmp.str[j], sep=":")
        }
      }
      ### Finish QLT name ###
                          
      # find the indices of the term(s) to be dropped
      # All terms contain label.term.drop will be dropped
      tmp.str <- strsplit(label.term.drop,":")[[1]]
      idx.term.drop <- 1:length(f.order)
      for(j in 1:length(tmp.str)) 
        idx.term.drop <- intersect(idx.term.drop,
                                   grep(tmp.str[j], f.label, ignore.case=TRUE))
      # the indices of term(s) to be kept
      idx.term.kept <- setdiff(1:length(f.order), idx.term.drop)
      
      #### regenerate a formula with the kept terms additive ###
      if(length(idx.term.kept) == 0) { # nothing left after drop label.term.drop
        msg <- paste("There will be nothing left if drop ", drop.term.name[i])
        stop(msg)
      }
      else {
        # All terms for idx.term.kept will be additive
        # Why it's so awkward? paste can't concatenate a list of strings?
        formula.new <- NULL
        for(j in 1:length(idx.term.kept)) {
          formula.new <- paste(formula.new, f.label[idx.term.kept[j]], sep="+")
        }
        formula.new <- as.formula(paste("y~", substr(formula.new, 2, 100000), sep=""))
      }
      ### Finish generating a new formula

      ### Start fitting model again
      # parse the input formula
      p.new <- parseformula(formula.new, dimnames(qtl$geno)[[2]], dimnames(covar)[[2]])
      n.gen.QC <- c(n.gen[p.new$idx.qtl]-1, rep(1, p.new$n.covar))

      # covarariate to be passed to C function
      covar.C <- NULL
      if(!is.null(p.new$idx.covar))
        covar.C <- as.matrix(covar[,p.new$idx.covar])
      
      # call C function to do the genome scan
      if(method == "imp") {
        z <- .C("R_fitqtl_imp",
                as.integer(n.ind), # number of individuals
                as.integer(p.new$n.qtl), # number of qtls
                as.integer(n.gen.QC), # number of genotypes QTLs and covarariates
                as.integer(n.draws), # number of draws
                as.integer(qtl$geno[,p.new$idx.qtl,]), # genotypes for selected marker
                as.integer(p.new$n.covar), # number of covarariate
                as.double(covar.C), # covarariate
                as.integer(p.new$formula.intmtx),  # formula matrix for interactive terms
                as.integer(p.new$n.int), # number of interactions in the formula
                as.double(pheno), # phenotype
                # return variables
                lod=as.double(0), # LOD score
                df=as.integer(0), # degree of freedom
                
                PACKAGE="qtl")
      }

      # record the result for dropping this term
      # df
      result[i,1] <- result.full[1,1] - z$df
      # LOD score
      result[i,3] <- result.full[1,4] - z$lod
      # % variance explained
      result[i,4] <- result.full[1,5] - 100*(1 - exp(-2*z$lod*log(10)/n.ind))
      # Type III SS for this term - computed from %var
      result[i,2] <- result.full[3,2] * result[i,4] / 100
      # P value (chi2)
      result[i,5] <- 1 - pchisq(2*log(10)*result[i,3], result[i,1])
      # P value (F)
      df0 <- length(pheno) - z$df - 1; df1 <- result.full[2,1];
      Rss0 <- result.full[2,2] + result[i,2];
      Rss1 <- result.full[2,2]
      Fstat <- ((Rss0-Rss1)/(df0-df1)) / (Rss1/df1)
      result[i,6] <- 1 - pf(Fstat, df0-df1, df1)

      rownames(result)[i] <- drop.term.name[i]
    } # finish dropping terms loop

    # assign output object
    output$result.drop <- result
    
  }  ## if(dropone)
      
#  else {
    # don't do drop one at at time
    # output the lod, pvar and df for this model
#    result <- matrix(rep(0,4),1,4)
#    result[1] <- z$lod
#    result[2] <- 100*(1 - exp(-2*z$lod*log(10)/n.ind))
#    result[3] <- z$df
#    result[4] <- 1 - pchisq(2*log(10)*z$lod, z$df)
#    rownames(result) <- "Full"
#    colnames(result) <- c("LOD", "%var", "df", "Pvalue")
#  }

  
  class(output) <- "fitqtl"
  attr(output, "method") <- method
  attr(output, "formula") <- formula
  attr(output, "type") <- qtl$type
  attr(output, "nind") <- length(pheno)
  output

}


#####################################################################
#
# parseformula
#
# Function to be called by fitqtl. It's used to
# parse the input formula
#
# This is the internal function and not supposed to be used by user
#
#####################################################################

parseformula <- function(formula, qtl.dimname, covar.dimname)
{
  # The terms for input formula
  f.formula <- terms(formula)
  order.term <- attr(f.formula, "order") # get the order of the terms
  idx.term <- which(order.term==1) # get the first order terms
  label.term <- attr(f.formula, "term.labels")[idx.term]
  formula.mtx <- attr(f.formula, "factors") # formula matrix
 
  idx.qtl <- NULL
  idx.covar <- NULL

  # loop thru all terms and find out how many QTLs and covarariates
  # are there in the formula. Construct idx.qtl and idx.covar at the same time
  for (i in 1:length(idx.term)) {
    # find out if there term is a QTL or a covarariate
    # ignore the case for QTLs, e.g., Q1 is equivalent to q1
    idx.tmp <- grep(paste(label.term[i],"$", sep=""),
                    qtl.dimname, ignore.case=TRUE)
    if( length(idx.tmp) )  # it's a QTL
      idx.qtl <- c(idx.qtl, idx.tmp)
    else if(label.term[i] %in% covar.dimname) # it's a covarariate
      idx.covar <- c(idx.covar, which(label.term[i]==covar.dimname))
    else
      stop(paste("Unrecognized term", label.term[i], "in formula"))
  }
  n.qtl <- length(idx.qtl) # number of QTLs in formula
  n.covar <- length(idx.covar) # number of covarariates in formula
  # now idx.qtl and idx.covar are the indices for genotype
  # and covarariate matrices according to input formula
 
  # loop thru all terms again and reorganize formula.mtx
  formula.idx <- NULL
  ii <- 1
  jj <- 1
  for (i in 1:length(idx.term)) {
    if(label.term[i] %in% qtl.dimname) {  # it's a QTL
      formula.idx <- c(formula.idx, ii)
      ii <- ii+1
    }
    else { # it's a covarariate
      formula.idx <- c(formula.idx, jj+n.qtl)
      jj <- jj+1
    }
  }

  # reorganize formula.mtx according to formula.idx
  # remove the first row (for y)
  formula.mtx <- formula.mtx[2:dim(formula.mtx)[1],]
  # rearrange the rows according to formula.idx if there's more than one row
  if(length(formula.idx) > 1)
    formula.mtx <- formula.mtx[order(formula.idx),]
  # take out only part of the matrix for interactions and pass to C function
  # all the input QTLs and covarariates for C function will be additive
  n.int <- length(order.term) - length(idx.term) # number of interactions
  if(n.int != 0)
    formula.intmtx <- formula.mtx[,(length(idx.term)+1):length(order.term)]
  else # no interaction terms
    formula.intmtx <- NULL
  

  # return object
  result <- NULL
  result$idx.qtl <- idx.qtl
  result$n.qtl <- n.qtl
  result$idx.covar <- idx.covar
  result$n.covar <- n.covar
  result$formula.intmtx <- formula.intmtx
  result$n.int <- n.int

  result

}


#####################################################################
#
# summary.fitqtl
#
#####################################################################
summary.fitqtl <- function(object, ...)
{
  # this is just an interface.
  class(object) <- "summary.fitqtl"
  object
}


#####################################################################
#
# print.summary.fitqtl
#
#####################################################################
print.summary.fitqtl <- function(x, ...)
{
  cat("\n")
  cat("\t\tSummary for fit QTL\n\n")
  cat( paste("Method is: ", attr(x, "method"), "\n") )
  cat( paste("Number of observations: ", attr(x, "nind"), "\n\n") )

  # print ANOVA table for full model
  cat("Full model result\n")
  cat("----------------------------------  \n")
  cat( paste("Model formula is: ", deparse(attr(x, "formula")), "\n\n") )
  print(x$result.full, quote=FALSE, na.print="")
  cat("\n\n")
  
  # print ANOVA table for dropping one at a time analysis (if any)
  if("result.drop" %in% names(x)) {
    cat("Drop one QTL at a time ANOVA table: \n")
    cat("----------------------------------  \n")
    # use print.coefmat instead of print.data.frame
    # make sure the last column is P value
    print.coefmat(x$result.drop, digits=4, cs.ind=1, P.values=TRUE, has.Pvalue=TRUE)
    cat("\n")
  }

}
