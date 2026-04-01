##--------------------------------------------
## required packages
message("Load packages")
suppressPackageStartupMessages({
    library(rmarkdown)
    library(knitr)
    library(devtools)
    library(yaml)
    library(BBmisc)
    library(GenomicAlignments)
    library(tidyr)
    library(data.table)
    library(dplyr)
    library(plotly)
    library(DelayedMatrixStats)
    library(FRASER)
    library(rhdf5)
})

contains_symlink <- function(paths) {
    symlinks <- Sys.readlink(list.files(paths, include.dirs=TRUE, recursive=TRUE, full.names=TRUE))
    return(length(symlinks[symlinks != ""]) > 0)
}

## Copy the savedObjects directory if it exists
saved_objects <- "./input/savedObjects"
if (dir.exists(saved_objects)) {
    if(contains_symlink(saved_objects)) {
        # Copy the savedObjects dir if it contains symlinks
        print("Copying savedObjects directory...")
        system(paste("cp", "-Lr", saved_objects, "savedObjects"))
    } else {
        # Symlink the savedObjects dir if it doesn't contains symlinks
        print("Symlinking savedObjects directory...")
        system(paste("ln", "-s", saved_objects, "savedObjects"))
    }
}

## Copy the cache directory if it exists
cache_dir <- "./input/cache"
if (dir.exists(cache_dir)) {
    if(contains_symlink(cache_dir)) {
        # Copy the cache dir if it contains symlinks
        print("Copying cache directory...")
        system(paste("cp", "-Lr", cache_dir, "cache"))
    } else {
        # Symlink the cache dir if it doesn't contains symlinks
        print("Symlinking cache directory...")
        system(paste("ln", "-s", cache_dir, "cache"))
    }
}

## helper functions
write_tsv <- function(x, file, row.names = FALSE, ...){
    write.table(x=x, file=file, quote=FALSE, sep='\t', row.names= row.names, ...)
}

extract_params <- function(params) {
    unlist(params)[1]
}

options("FRASER.maxSamplesNoHDF5"=0)
options("FRASER.maxJunctionsNoHDF5"=-1)

h5disableFileLocking()

# set psiTypes to run based on fraser version
configure_fraser <- function(fraser_version) {
    if(fraser_version == "FRASER2"){
        pseudocount(0.1)
        psiTypes <- c("jaccard")
        psiTypesNotUsed <- c("psi5", "psi3", "theta")
    } else{
        pseudocount(1)
        psiTypes <- c("psi5", "psi3", "theta")
        psiTypesNotUsed <- c("jaccard")
    }
}
