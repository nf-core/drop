#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/mae-pipeline/MAE/gene_name_mapping.R

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')


suppressPackageStartupMessages({
    library(rtracklayer)
    library(data.table)
    library(magrittr)
    library(tidyr)
})

gtf_dt <- import("$gtf") %>% as.data.table
if (!"gene_name" %in% colnames(gtf_dt)) {
    gtf_dt[gene_name := gene_id]
}
if('gene_biotype' %in% colnames(gtf_dt))
    setnames(gtf_dt, 'gene_biotype', 'gene_type')
gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type)]

# make gene_names unique
gtf_dt[, N := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_orig := gene_name]
gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]

fwrite(gtf_dt, paste0("gene_name_mapping_", prefix, ".tsv"), na = NA)


## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
        paste('    bioconductor-rtracklayer:', as.character(packageVersion('rtracklayer')))
    ),
'versions.yml')
