
Sample data for the R/qtl package

These files contain sample data in several formats, so that the user
may better understand how data may be formatted for import via the
read.cross function.  These are the same as the "listeria" data
set included with the R/qtl package.  

Note: Replace the "..." in the directory string to the appropriate
location of the sampledata directory (for example,
"/usr/local/lib/R/library/qtl" or "c:/R/rw1051/library/qtl").


1. "csv" format

    File:

        listeria.csv

    Data import:

        listeria.a <- read.cross("csv", ".../sampledata", "listeria.csv")

2. "mm" (mapmaker) format

   Files:

       listeria.raw    "raw" file (with genotype and phenotype data)
       listeria.map    Genetic map information (markers must be in
                       order; map positions are not required)
       listeria.maps   Genetic map information, as produced by 
                       Mapmaker/exp

   Data import:

       listeria.b <- read.cross("mm", ".../sampledata",
                                "listeria.raw",,"listeria.map")

       listeria.bb <- read.cross("mm", ".../sampledata",
                                 "listeria.raw",,"listeria.maps")


3. "qtx" (Mapmanager QTX) format

   File:

       listeria.qtx

   Data import:

       listeria.c <- read.cross("qtx", ".../sampledata", "listeria.qtx")


4. "qtlcart" (QTL Cartographer) format

   Files:

        listeria_qc.cro    Genotype/phenotype data
        listeria_qc.map    Genetic map information

   Data import:

        listeria.d <- read.cross("qtlcart", ".../sampledata", 
                                 "listeria_qc.cro", "listeria_qc.map")


5. "karl" format

    Files:

        gen.txt    Genotype data
        phe.txt    Phenotype data
        map.txt    Genetic map information

    Data import:

        listeria.e <- read.cross("karl", ".../sampledata",
                                 genfile="gen.txt", phefile="phe.txt",
                                 mapfile="map.txt")

    or just the following:  

        listeria.e <- read.cross("karl", ".../sampledata")
