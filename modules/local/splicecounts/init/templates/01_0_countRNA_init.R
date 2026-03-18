#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_0_SPLICECOUNTS_INIT.R

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

source("$config", echo=FALSE)
if("$fraser_version" == "FRASER2"){
    pseudocount(0.1)
    psiTypes <- c("jaccard")
    psiTypesNotUsed <- c("psi5", "psi3", "theta")
} else{
    pseudocount(1)
    psiTypes <- c("psi5", "psi3", "theta")
    psiTypesNotUsed <- c("jaccard")
}

dataset    <- "$drop_group"
colDataFile <- "$col_data"
workingDir <- "./"

# Create initial FRASER object
col_data <- fread(colDataFile)

col_data\$strand <- 0L

fds <- FraserDataSet(colData = col_data,
                    workingDir = workingDir,
                    name       = paste0("raw-local-", dataset))

# Add paired end and strand specificity to the fds
pairedEnd(fds) <- colData(fds)\$pairedEnd
strandSpecific(fds) <- colData(fds)\$strand

# Save initial FRASER dataset
fds <- saveFraserDataSet(fds)

message(date(), ": FRASER object initialized for ", dataset)

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
        paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
    ),
'versions.yml')
