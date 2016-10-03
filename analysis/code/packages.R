## Add any new packages to the `pkgs' vector in alphabetical order.
pkgs <- c("knitr", "secr", "testthat", "xtable")
for (i in pkgs){
    if (!require(i, character.only = TRUE)){
        install.packages(i, repos = "http://star-www.st-andrews.ac.uk/cran/")
    }
}
