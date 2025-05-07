#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/template/Scripts/Pipeline/preprocessGeneAnnotation.R

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(data.table)
    library(rtracklayer)
    library(magrittr)
})

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

## Create txdb
txdb <- makeTxDbFromGFF("$gtf")
txdb <- keepStandardChromosomes(txdb)

tmpFile <- tempfile()
saveDb(txdb, tmpFile)
R.utils::copyFile(tmpFile, paste(prefix, ".db", sep=""), overwrite=TRUE)

# save count ranges
count_ranges <- exonsBy(txdb, by = "gene")
saveRDS(count_ranges, paste(prefix, ".count_ranges.Rds", sep=""))

# Get a gene annotation table
gtf_dt <- import("$gtf") %>% as.data.table
if (!"gene_name" %in% colnames(gtf_dt)) {
    gtf_dt[, gene_name := gene_id]
}
gtf_dt <- gtf_dt[type == 'gene']

if('gene_biotype' %in% colnames(gtf_dt))
    setnames(gtf_dt, 'gene_biotype', 'gene_type')

# Subset to the following columns only
columns <- c('seqnames', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type', 'gene_status')
columns <- intersect(columns, colnames(gtf_dt))
gtf_dt <- gtf_dt[, columns, with = FALSE]

# make gene_names unique
gtf_dt[, N := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_orig := gene_name]
gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]

fwrite(gtf_dt, paste(prefix, ".gene_name_mapping.tsv", sep=""), na = NA)

# Write the versions file

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    bioconductor-genomicscores:', as.character(packageVersion('GenomicFeatures'))),
        paste('    bioconductor-rtracklayer:', as.character(packageVersion('rtracklayer'))),
        paste('    bioconductor-txdbmaker:', as.character(packageVersion('txdbmaker')))
    ),
'versions.yml')
