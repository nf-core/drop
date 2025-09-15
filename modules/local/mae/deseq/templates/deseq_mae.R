#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/mae-pipeline/MAE/deseq_mae.R

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

suppressPackageStartupMessages({
    library(stringr)
    library(tMAE)
})

message("Started with deseq")

# Read mae counts
mae_counts <- fread("$counts", fill=TRUE)
mae_counts <- mae_counts[contig != '']
mae_counts[, position := as.numeric(position)]

# Sort by chr
mae_counts <- mae_counts[!grep("Default|opcode", contig)]
mae_counts[,contig := factor(contig,
                levels = unique(str_sort(mae_counts\$contig, numeric = TRUE)))]


print("Running DESeq...")
# Function from tMAE pkg
rmae <- DESeq4MAE(mae_counts) ## negative binomial test for allelic counts

### Add AF information from gnomAD
if (${add_af ? 'TRUE' : 'FALSE'} == TRUE) {
    print("Adding gnomAD allele frequencies...")

    # obtain the assembly from the config
    genome_assembly <- "${assembly}"
    rmae <- add_gnomAD_AF(rmae, genome_assembly = genome_assembly,
        max_af_cutoff = $max_af, populations = c("AF", "AF_afr", "AF_amr", "AF_eas", "AF_nfe"))
} else {
    rmae[, rare := NA]
}

saveRDS(rmae, paste0(prefix, "_res.Rds"))

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-r.utils:', as.character(packageVersion('R.utils'))),
        paste('    r-string:', as.character(packageVersion('stringr'))),
        paste('    r-tmae:', as.character(packageVersion('tMAE')))
    ),
'versions.yml')
