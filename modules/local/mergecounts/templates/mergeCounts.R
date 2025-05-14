#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-expression-pipeline/Counting/mergeCounts.R

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(BiocParallel)
    library(SummarizedExperiment)
})

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

register(MulticoreParam($task.cpus))
count_ranges <- readRDS("$count_ranges")

# Read counts
counts_list <- bplapply(list(${counts.collect { file -> "\"$file\""}.join(', ')}), function(f){
    if(grepl('Rds\$', f))
        assay(readRDS(f))
    else {
        ex_counts <- as.matrix(fread(f), rownames = "geneID")
        print(head(ex_counts))
        stopifnot(! list(${exclude_ids.collect { id -> "\"$id\""}.join(', ')}) %in% names(ex_counts))
        subset(ex_counts, select = ${exclude_ids.collect { id -> "\"$id\""}.join(', ')})
    }
})
message(paste("read", length(counts_list), 'files'))

# check rownames and proceed only if they are the same
row_names_objects <- lapply(counts_list, rownames)
if( length(unique(row_names_objects)) > 1 ){
    stop('The rows (genes) of the count matrices to be merged are not the same.')
}

# merge counts
merged_assays <- do.call(cbind, counts_list)
total_counts <- SummarizedExperiment(assays=list(counts=merged_assays))
colnames(total_counts) <- gsub('.bam', '', colnames(total_counts))

# assign ranges
rowRanges(total_counts) <- count_ranges

# Add sample annotation data (colData)
sample_anno <- fread("$samplesheet",
                    colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
sample_anno <- sample_anno[, .SD[1], by = RNA_ID]
col_data <- data.table(RNA_ID = as.character(colnames(total_counts)))
col_data <- left_join(col_data, sample_anno, by = "RNA_ID")
rownames(col_data) <- col_data\$RNA_ID
colData(total_counts) <- as(col_data, "DataFrame")
rownames(colData(total_counts)) <- colData(total_counts)\$RNA_ID

# save in RDS format
saveRDS(total_counts, paste(prefix, ".Rds", sep=""))

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    r-data.table:', as.character(packageVersion('dplyr'))),
        paste('    bioconductor-biocparallel:', as.character(packageVersion('BiocParallel'))),
        paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment')))
    ),
'versions.yml')
