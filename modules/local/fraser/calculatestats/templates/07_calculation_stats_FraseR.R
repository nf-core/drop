#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/07_calculation_stats_FraseR.R

source(Sys.which("aberrant_splicing_config.R"), echo=FALSE)
if("$fraser_version" == "FRASER2"){
    pseudocount(0.1)
    psiTypes <- c("jaccard")
    psiTypesNotUsed <- c("psi5", "psi3", "theta")
} else{
    pseudocount(1)
    psiTypes <- c("psi5", "psi3", "theta")
    psiTypesNotUsed <- c("jaccard")
}

annotation <- "$annotation_id"
dataset    <- "$drop_group"
fdsFile    <- "$fds"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# read in subsets from sample anno if present (returns NULL if not present)
source(Sys.which("parse_subsets_for_FDR.R"))
fraser_sample_ids <- list(${sample_ids.collect { "\"$it\"" }.join(', ')})
subsets <- parse_subsets_for_FDR("$genes_to_test",
                                sampleIDs=fraser_sample_ids)

# Load FRASER data
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# Calculate stats
for (type in psiTypes) {
    # Pvalues
    fds <- calculatePvalues(fds, type=type)
    # Adjust Pvalues
    fds <- calculatePadjValues(fds, type=type, subsets=subsets)
}

fds <- saveFraserDataSet(fds)

# remove .h5 files from previous runs with other FRASER version
fdsDir <- "$fds/${drop_group}--${annotation_id}/"
pvalFiles <- grep("p(.*)BetaBinomial_(.*).h5",
                list.files(fdsDir),
                value=TRUE)
for(type in psiTypesNotUsed){
    pvalFilesType <- grep(type, pvalFiles, value=TRUE)
    for(pFile in pvalFilesType){
        unlink(file.path(fdsDir, pFile))
    }
}

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
