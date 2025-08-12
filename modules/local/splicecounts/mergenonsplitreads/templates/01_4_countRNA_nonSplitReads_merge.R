#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_4_countRNA_nonSplitReads_merge.R

source("$config", echo=FALSE)
configure_fraser("$fraser_version")

dataset    <- "$drop_group"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Move non split counts to the correct cache directory
system("mv cache/nonSplicedCounts/raw-local cache/nonSplicedCounts/raw-local-${drop_group}")

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS("$non_split_counts_granges")

# Remove the cache file of merged counts from previous runs
nonSplitCountsDir <- file.path(workingDir, "savedObjects",
                            paste0("raw-local-", dataset), 'nonSplitCounts')
if(dir.exists(nonSplitCountsDir)){
    unlink(nonSplitCountsDir, recursive = TRUE)
}

# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds,
                                                    splitCountRanges=splitCounts_gRanges,
                                                    minAnchor=5,
                                                    recount=FALSE,
                                                    longRead=${long_read ? "TRUE" : "FALSE"},)

# remove cache dir as its not needed later
if(dir.exists(nonSplitCountsDir)){
    unlink(nonSplitCountsDir, recursive = TRUE)
}

message(date(), ":", dataset, " nonSplit counts done")

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
