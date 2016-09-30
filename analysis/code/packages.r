## Add any new packages to the `pkgs' vector.
pkgs <- "secr"
for (i in pkgs){
    if (!require(i, character.only = TRUE)){
        install.packages(i, repos = "http://star-www.st-andrews.ac.uk/cran/")
    }
}
