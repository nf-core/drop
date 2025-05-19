#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(magrittr)
    library(tools)
    library(yaml)
})

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

ods <- readRDS("$ods_fitted")
implementation <- "$implementation"
register(MulticoreParam($task.cpus))

# read in gene subsets from sample anno if present (returns NULL if not present)
source("$parse_subsets_for_FDR")
outrider_sample_ids <- list(${ids.collect { return "\"$it\"" }.join(', ')})
subsets <- parse_subsets_for_FDR("$genes_to_test", sampleIDs=outrider_sample_ids)
subsets_final <- convert_to_geneIDs(subsets, "$gene_name_mapping")

# P value calculation
message(date(), ": P-value calculation ...")
ods <- computePvalues(ods, subsets=subsets_final)
message(date(), ": Zscore calculation ...")
ods <- computeZscores(ods, peerResiduals=grepl('^peer\$', implementation))

# save ods with pvalues
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
        paste('    r-tools:', as.character(packageVersion('tools'))),
        paste('    r-yaml:', as.character(packageVersion('yaml'))),
        paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment')))
    ),
'versions.yml')
