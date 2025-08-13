#!/usr/bin/env Rscript
#https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-expression-pipeline/OUTRIDER/runOutrider.R

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(magrittr)
    library(tools)
})

ods <- readRDS("$ods")
implementation <- "$implementation"
mp <- $ae_max_tested_dimension_proportion
register(MulticoreParam(${task.cpus}))

## subset filtered
ods <- ods[mcols(ods)\$passedFilter,]

# add gene ranges to rowData
gr <- unlist(endoapply(rowRanges(ods), range))
if(length(gr) > 0){
    rd <- rowData(ods)
    rowRanges(ods) <- gr
    rowData(ods) <- rd
}

ods <- estimateSizeFactors(ods)

## find optimal encoding dimension
a <- 5
b <- min(ncol(ods), nrow(ods)) / mp   # N/3

maxSteps <- 15
if(mp < 4){
    maxSteps <- 20
}

Nsteps <- min(maxSteps, b)   # Do at most 20 steps or N/3
# Do unique in case 2 were repeated
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
ods <- findEncodingDim(ods, params = pars_q, implementation = implementation)
opt_q <- getBestQ(ods)

## fit OUTRIDER
# ods <- OUTRIDER(ods, implementation = implementation)
message(date(), ": SizeFactor estimation ...")
ods <- estimateSizeFactors(ods)
message(date(), ": Controlling for confounders ...")
implementation <- tolower(implementation)
ods <- controlForConfounders(ods, q=opt_q, implementation=implementation)
if (implementation == "peer" || implementation == "pca"){
    message(date(), ": Fitting the data ...")
    ods <- fit(ods)
}
message("outrider fitting finished")

saveRDS(ods, paste(prefix, ".Rds", sep=""))

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
        paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
        paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
        paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
    ),
'versions.yml')
