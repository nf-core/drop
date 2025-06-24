#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/04_fit_hyperparameters_FraseR.R

source("$config", echo=FALSE)
configure_fraser("$fraser_version")

rseed <- ${random_seed instanceof Integer ? random_seed : random_seed ? 'TRUE' : 'FALSE'}
if (!isFALSE(rseed)) {
    if(isTRUE(rseed)){
        set.seed(42)
    } else if (is.numeric(rseed)){
        set.seed(as.integer(rseed))
    }
}

#+ input
dataset    <- "$drop_group"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Load PSI data
fds <- loadFraserDataSet(dir=workingDir, name=dataset)
fitMetrics(fds) <- psiTypes

# Run hyper parameter optimization
implementation <- "$implementation"
mp <- $max_tested_dimension_proportion

# Get range for latent space dimension
a <- 2
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

maxSteps <- 12
if(mp < 6){
    maxSteps <- 15
}

Nsteps <- min(maxSteps, b)
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique

for(type in psiTypes){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type,
                            implementation=implementation,
                            q_param=pars_q,
                            plot = FALSE)
    fds <- saveFraserDataSet(fds)
}
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
