
README file for the qtl package

copyright (c) 2001, Karl W Broman, Johns Hopkins University
Licensed under the GNU General Public License version 2 (June, 1991)

This is a quick document to explain the installation and use
of the "qtl" package for R


Installation (Windows)

  1. Unzip the "qtl.zip" file into the directory $RHOME\library
     ($RHOME is something like c:\Program Files\R\rw1030)
     Note that this should create a directory $RHOME\library\qtl
     containing the R source code and the compiled dll

  2. Start Rgui 

  3. Type "link.html.help()" to get the help files for the qtl 
     package added to the help indices

  4. Note that the source code is in the file "qtl_*.tar.gz"


Installation (Unix)

  1. Go into the directory containing "qtl_*.tar.gz"

  2. Type "R INSTALL qtl"


A brief example of using the package 
(See the help file "A starting point" for further information.)

library()                  # List available packages
library(help=qtl)          # List functions and data within "qtl"
library(qtl)               # Load the qtl package

data(listeria)             # Make the dataset "listeria" available

# different ways to get help
?sim.cross
help(sim.cross)
help.start()

plot.missing(listeria)     # plot genotype matrix
plot.map(listeria)         # plot genetic map
plot(listeria)             # both of these, plus phenotype histograms 

# remove markers with no genotype data
listeria <- remove.nullmarkers(listeria)

# estimate genetic map
newmap <- est.map(listeria,error.prob=0,print.rf=TRUE)    
plot.map(listeria,newmap)       # plot original map vs new map

# calc cond'l genotype probabilities
listeria <- calc.genoprob(listeria,step=1,off.end=0)  

# identify most likely underlying genotypes
listeria <- argmax.geno(listeria)

# identify genotypes possibly in error
listeria <- calc.errorlod(listeria)
plot.errorlod(listeria)
top.errorlod(listeria)

# estimate pairwise recombination fractions
listeria <- est.rf(listeria)
plot.rf(listeria)
plot.rf(listeria,c(1,5,13))

# take log phenotype and then do genome scan
# (for these data, this is not really appropriate)
listeria$pheno[,1] <- log(listeria$pheno[,1])
output1 <- scanone(listeria,method="anova")
output2 <- scanone(listeria,method="im")
plot(output1)
plot(output1,c(5,13))
plot(output1,output2,c(5,13))

