#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/mae-pipeline/QC/deseq_qc.R

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(tMAE)
})

# Read MA counts for qc
qc_counts <- fread("$counts", fill=TRUE)
# qc_counts <- qc_counts[!is.na(position)]

# Run DESeq
rmae <- DESeq4MAE(qc_counts, minCoverage = 10)
rmae[, RNA_GT := '0/1']
rmae[altRatio < .2, RNA_GT := '0/0']
rmae[altRatio > .8, RNA_GT := '1/1']
rmae[, position := as.numeric(position)]

# Convert to granges
qc_gr <- GRanges(seqnames = rmae\$contig,
                ranges = IRanges(start = rmae\$position, end = rmae\$position),
                strand = '*')
mcols(qc_gr) = DataFrame(RNA_GT = rmae\$RNA_GT)

saveRDS(qc_gr, paste0(prefix, ".Rds"))

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-tmae:', as.character(packageVersion('tMAE'))),
        paste('    bioconductor-genomicranges:', as.character(packageVersion('GenomicRanges')))
    ),
'versions.yml')
