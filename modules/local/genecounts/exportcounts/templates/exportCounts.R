#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-expression-pipeline/Counting/exportCounts.R

suppressPackageStartupMessages({
    library(data.table)
    library(SummarizedExperiment)
})

total_counts <- readRDS("${counts}")

# save in exportable format
fwrite(as.data.table(assay(total_counts), keep.rownames = 'geneID'),
    file = "geneCounts.tsv.gz",
    quote = FALSE, row.names = FALSE, sep = '\t', compress = 'gzip')

# Write the versions file
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
        paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
    ),
'versions.yml')
