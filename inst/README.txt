
README file for the qtl package
----------------------------------------------------------------------
This explains the installation of the qtl package for R.  If you have
any problems with installation, send an email to Karl Broman
<kbroman@jhsph.edu>.
----------------------------------------------------------------------

Obtaining R

  You can download R from the Comprehensive R Archive Network (CRAN).
  Visit http://cran.r-project.org or a local mirror (for example,
  http://cran.us.r-project.org).  Source code is available for UNIX,
  and binaries are available for Windows, MacOS, and many versions of
  Linux.  

  On Windows, we recommend installing R in "c:\R" rather than
  "c:\Program Files\R".  Why didn't Microsoft use "Programs" rather
  than "Program files"?


Obtaining R/qtl

  You can obtain the latest version of R/qtl from 

      http://biosun01.biostat.jhsph.edu/~kbroman/qtl

  Copies of R/qtl will also be placed on CRAN, but the version at the 
  above site will be update more frequently.


Installation of R/qtl (Windows)

  1. Unzip the "qtl.zip" file into the directory $RHOME\library (where
     $RHOME is something like c:\R\rw1030).  Note that this should
     create a directory $RHOME\library\qtl containing the R source
     code and the compiled dll.  

  2. Start Rgui. 

  3. Type "link.html.help()" to get the help files for the qtl 
     package added to the help indices.

  4. Note that the source code is in the file "qtl_*.tar.gz", which
     you may download separately.


Installation of R/qtl (Unix)

  1. Go into the directory containing "qtl_*.tar.gz".

  2. Type "R INSTALL qtl" to have the package installed in the
     standard location (/usr/local/lib/R/library).  You'll probably
     need to be superuser.

     To install the package locally, type 

         R INSTALL --library=$LOCALLIB qtl_*.tar.gz

     (where $LOCALLIB is something like "/home/auser/Rlibs").  You'll
     then need to create a file ~/.Renviron containing the line 

         R_LIBS=/home/auser/Rlibs

     so that R will know to search for packages in that directory.


Getting started

  Once you start R, you'll need to type "library(qtl)" to load the
  package.  You can create a file "~/.Rprofile" or "c:\.Rprofile"
  containing R code to be run whenever you start R.  If you use the
  R/qtl package regularly, you should place the line "library(qtl)" in
  such a file.

  Efficient use of the R/qtl package requires considerable knowledge
  of the R language.  Learning R may require a formidable investment
  of time, but it will definitely be worth the effort.  Numerous free
  documents on getting started with R are available on CRAN
  (http://cran.r-project.org).  In additional, several books are
  available on S/S-PLUS, which is very similar to R.  For example, see
  WN Venables, BD Ripley (1999) Modern Applied Statistics with S-PLUS,
  3rd edition. Springer.

  To get started with R/qtl, you might first peruse the documentation
  that is bundled with it.  Type help.start() to start the html
  version of the R help pages.  Then click "Packages" -> "qtl".  

  In Windows, you may gain access to the help documents by clicking
  "Help" in the menu bar and then "R language (html)".  Windows users
  may wish to place the line "options(htmlhelp=TRUE)" in the
  file "c:\.Rprofile".  

  The help file titled "A starting point" gives a brief walk-through
  of an example analysis, and so is a good place to start.  You may
  also view this help file by typing 

      ?"A starting point"

  from the command line in R.  

  A tutorial on R/qtl (as a PDF document) is also available.  It
  briefly describes the aims of the R/qtl package, lists the available
  functions grouped in categories, and provides a few extended
  examples.  The tutorial is bundled with R/qtl, as "rqtltour.pdf" and
  is also available from the R/qtl website: 

      http://biosun01.biostat.jhsph.edu/~kbroman/qtl


Questions/comments/concerns

  If you have any questions, suggestions, problems or complaints
  regarding R/qtl, please email Karl Broman <kbroman@jhsph.edu>.
  
----------------------------------------------------------------------
end of README.txt
