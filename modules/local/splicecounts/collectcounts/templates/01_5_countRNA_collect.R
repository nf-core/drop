#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_5_SPLICECOUNTS_COLLECTCOUNTS.R

source("$config", echo=FALSE)
configure_fraser("$fraser_version")

dataset    <- "$drop_group"
workingDir <- "./"
saveDir    <- "savedObjects/raw-local-$drop_group"

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))
splitCounts_gRanges <- readRDS("$splitcounts_granges")
spliceSiteCoords <- readRDS("$splicesites_splitcounts")

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts_h5 <- HDF5Array::HDF5Array(file.path(saveDir, "rawCountsJ.h5"), "rawCountsJ")
splitCounts_se <- SummarizedExperiment(
    colData = colData(fds),
    rowRanges = splitCounts_gRanges,
    assays = list(rawCountsJ=splitCounts_h5)
)


nonSplitCounts_h5 <- HDF5Array::HDF5Array(file.path(saveDir, "rawCountsSS.h5"), "rawCountsSS")
nonSplitCounts_se <- SummarizedExperiment(
    colData = colData(fds),
    rowRanges = spliceSiteCoords,
    assays = list(rawCountsSS=nonSplitCounts_h5)
)

# Add Counts to FRASER object
fds <- addCountsToFraserDataSet(fds=fds, splitCounts=splitCounts_se,
                                nonSplitCounts=nonSplitCounts_se)

# Save final FRASER object
fds <- saveFraserDataSet(fds)

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-rmarkdown:', as.character(packageVersion('rmarkdown'))),
        paste('    r-knitr:', as.character(packageVersion('knitr'))),
        paste('    r-devtools:', as.character(packageVersion('devtools'))),
        paste('    r-yaml:', as.character(packageVersion('yaml'))),
        paste('    r-bbmisc:', as.character(packageVersion('BBmisc'))),
        paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
        paste('    r-plotly:', as.character(packageVersion('plotly'))),
        paste('    r-rhdf5:', as.character(packageVersion('rhdf5'))),
        paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
        paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
        paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
        paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
    ),
'versions.yml')
