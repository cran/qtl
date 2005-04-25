
README file for the qtl package
----------------------------------------------------------------------
This explains the installation of the qtl package for R.  If you have
any problems with installation, send an email to Karl Broman
<kbroman@jhsph.edu>.
----------------------------------------------------------------------

OBTAINING R

  You can download R from the Comprehensive R Archive Network (CRAN).
  Visit http://cran.r-project.org or a local mirror (for example,
  http://cran.us.r-project.org).  Source code is available for Unix,
  and binaries are available for Windows, MacOS, and many versions of
  Linux.  


OBTAINING R/QTL

  You can obtain the latest version of R/qtl from 

      http://www.biostat.jhsph.edu/~kbroman/qtl

  Copies of R/qtl will also be placed on CRAN (cran.r-project.org),
  but the version at the above site may be updated more frequently.
  Binaries are available for Windows and MacOS; source code is
  available for Unix.


INSTALLATION OF R AND R/QTL (Windows)

  1. The Windows version of R is distributed as a single file,
     with a name something like rw2010.exe.  Install R by executing
     this file.  We recommend installing R in "c:\R" rather than
     "c:\Program Files\R".  Why didn't Microsoft use "Programs" rather
     than "Program files"? 

  2. To install R/qtl, download the file "qtl.99-22.zip" (or the
     equivalent) and then do one of the following:

     a. Start R.  Select (on the menu bar) "Packages" and then
        "Install package from local zip file...".  Find the file
        "qtl_0.99-22.zip" on your hard drive, and click "Open".

     b. Unzip the "qtl_0.99-22.zip" file into the directory
        $RHOME\library (where $RHOME is something like c:\R\rw1081).
        Note that this should create a directory $RHOME\library\qtl
        containing the R source code and the compiled dll.

        Start R and type "link.html.help()" to get the help files for
        the qtl package added to the help indices.


INSTALLATION OF R AND R/QTL (MacOS version 10.2.x and above)

  1. Download the file RAqua.dmg and double-click it to mount a
     "drive" with a name something like "R-2.0.1.dmg".  Follow the
     instructions in the file "ReadMe.txt".

  2. To install R/qtl:

     a. Download the compiled version of R/qtl for Mac OS X, a file
        like "qtl_0.99-22.tgz".

     b. Start R by double-clicking "R" in your Applications folder.

     c. From the menu bar, click "Packages & Data" -> "Package
        Installer".  

	Select "Local Binary Package" from the drop-down menu at the
        top of the window that comes up.

	Click "Install/Update" at the bottom of the window.

	Find the package on your harddrive and click "Open".

	Close the window.


INSTALLATION OF R/QTL (Unix)

  1. We'll assume that R has already been installed. 

  2. Go into the directory containing the file "qtl_0.99-22.tar.gz".

  3. Do one of the following:

     a. To install R/qtl in the standard location
        (/usr/local/lib/R/library), type 

            R CMD INSTALL qtl_*.tar.gz

        You'll probably need to be superuser.

     b. To install the package locally, type 

            R CMD INSTALL --library=/home/auser/Rlibs qtl_*.tar.gz

        (where "/home/auser/Rlibs" should be replaced with the
        appropriate directory).  

        Create a file ~/.Renviron containing the line

            R_LIBS=/home/auser/Rlibs

        so that R will know to search for packages in that directory.


GETTING STARTED

  Once you start R, you'll need to type "library(qtl)" to load the
  package.  You can create a file "~/.Rprofile" (Unix or MacOS) or
  "c:\.Rprofile" (Windows) containing R code to be run whenever you
  start R.  If you use the R/qtl package regularly, you should place
  the line "library(qtl)" in such a file.

  Efficient use of the R/qtl package requires considerable knowledge
  of the R language.  Learning R may require a formidable investment
  of time, but it will definitely be worth the effort.  Numerous free
  documents on getting started with R are available on CRAN
  (http://cran.r-project.org).  In addition, several books are
  available.  For example, see WN Venables, BD Ripley (2002) Modern
  Applied Statistics with S, 4th edition. Springer.

  To get started with R/qtl, you might first peruse the documentation
  that is bundled with it.  Type help.start() to start the html
  version of the R help pages.  Then click "Packages" -> "qtl".  

  In Windows or MacOS, you may gain access to the help documents by
  clicking "Help" in the menu bar and then "R language (html)".  If
  you include "options(htmlhelp=TRUE)" in your .Rprofile file, use of
  the html version of the help pages will be automatic.

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

      http://www.biostat.jhsph.edu/~kbroman/qtl


CITING R/QTL

  To cite R/qtl in publications, use

      Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping
      in experimental crosses.  Bioinformatics 19:889-890


QUESTIONS/COMMENTS/CONCERNS

  If you have any questions, suggestions, problems or complaints
  regarding R/qtl, please email Karl Broman <kbroman@jhsph.edu>.
  
----------------------------------------------------------------------
end of README.txt
