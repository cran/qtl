
Sample data for the R/qtl package

These files contain sample data in several formats, so that the user
may better understand how data may be formatted for import via the
read.cross function.  These are the same as the "listeria" data
set included with the R/qtl package.


1. "csv" format

    File:

        listeria.csv

    Data import:

        listeria.a <- read.cross("csv", dir="../sampledata",
                                 file="listeria.csv")

2. "mm" (mapmaker) format

   Files:

       listeria.raw    "raw" file (with genotype and phenotype data)
       listeria.map    Genetic map information (markers must be in
                           order; map positions are not required

   Data import:

       listeria.b <- read.cross("mm", dir="../sampledata",
                                rawfile="listeria.raw",
                                mapfile="listeria.map")


3. "karl" format

    Files:

        gen.txt    Genotype data
        phe.txt    Phenotype data
        map.txt    Genetic map information

    Data import:

        listeria.c <- read.cross("karl", dir=".../sampledata",
                                 genfile="gen.txt", phefile="phe.txt",
                                 mapfile="map.txt")
