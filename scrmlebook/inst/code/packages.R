## Add any new packages to the `pkgs' vector in alphabetical order.
pkgs <- c("gdistance", "igraph", "fields", "knitr", "maptools", "raster", "secr", "sp", "spatstat", "testthat", "xtable")
args <- commandArgs(trailingOnly = TRUE)
upgrade <- "-u" %in% args
for (i in pkgs){
    if (!require(i, character.only = TRUE) | upgrade){
        install.packages(i, repos = "http://star-www.st-andrews.ac.uk/cran/")
    }
}
