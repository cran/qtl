######################################################################
#
# TestIO/input.R
#
# copyright (c) 2002, Karl W Broman, Johns Hopkins University
# last modified Feb, 2002
# first written Feb, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# This file contains code for testing the cross IO in R/qtl.
#
# Needed input files:
#
#    gen.txt, map.txt, phe.txt    [Karl's format]
#    listeria.raw, listeria.map   [mapmaker format]
#    listeria.raw, listeria2.map  [mapmaker format; no marker pos]
#    listeria.csv                 [csv format] 
#    listeria2.csv                [csv format; no marker pos] 
#
######################################################################

library(qtl)

##############################
# Reading
##############################
# Read Karl's format
karl <- read.cross("karl", genfile="gen.txt",
                   mapfile="map.txt", phefile="phe.txt")
karl
karl2 <- read.cross("karl")
karl2

# Read CSV format
csv <- read.cross("csv", , "listeria.csv")
csv
csv2 <- read.cross("csv", , "listeria2.csv")
csv2

# Read mapmaker format
mm <- read.cross("mm", "", "listeria.raw", "listeria.map")
mm
mm2 <- read.cross("mm", "", "listeria.raw", "listeria2.map")
mm2

##############################
# Writing
##############################
# Write in CSV format
write.cross(karl, "csv", filestem="junk1")
karl3 <- read.cross("csv", "", "junk1.csv")
comparecrosses(karl, karl3)

# Write in mapmaker format
write.cross(karl, "mm", filestem="junk2")
