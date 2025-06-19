#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_3_countRNA_nonSplitReads_samplewise.R

source("$config", echo=FALSE)

dataset    <- "$drop_group"
workingDir <- "./"

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Get sample id from wildcard
sample_id <- "$sample_id"


# Read splice site coordinates from RDS
spliceSiteCoords <- readRDS("${cache}/raw-local-${drop_group}/spliceSites_splitCounts.rds")

# Count nonSplitReads for given sample id
sample_result <- countNonSplicedReads(sample_id,
                                    splitCountRanges = NULL,
                                    fds = fds,
                                    NcpuPerSample = $task.cpus,
                                    minAnchor=5,
                                    recount=${recount ? "TRUE" : "FALSE"},
                                    spliceSiteCoords=spliceSiteCoords,
                                    longRead=${long_read ? "TRUE" : "FALSE"})

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (non split counts) = ", length(sample_result))

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
