#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/05_fit_autoencoder_FraseR.R

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
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

fds <- loadFraserDataSet(dir=workingDir, name=dataset)

# Fit autoencoder
# run it for every type
implementation <- "$implementation"

for(type in psiTypes){
    currentType(fds) <- type
    q <- bestQ(fds, type)
    verbose(fds) <- 3   # Add verbosity to the FRASER object
    fds <- fit(fds, q=q, type=type, iterations=15, implementation=implementation)
    fds <- saveFraserDataSet(fds)
}

# remove .h5 files from previous runs with other FRASER version
fdsDir <- "$fds/$drop_group/"
for(type in psiTypesNotUsed){
    predMeansFile <- file.path(fdsDir, paste0("predictedMeans_", type, ".h5"))
    if(file.exists(predMeansFile)){
        unlink(predMeansFile)
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
