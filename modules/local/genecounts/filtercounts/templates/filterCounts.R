#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-expression-pipeline/Counting/filterCounts.R

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

counts <- readRDS("$counts")
ods <- OutriderDataSet(counts)
txdb <- loadDb("$txdb")

# filter not expressed genes
fpkmCutoff <- $fpkmCutoff
ods <- filterExpression(ods, gtfFile=txdb, filter=FALSE,
                        fpkmCutoff=fpkmCutoff, addExpressedGenes=TRUE)

# add column for genes with at least 1 gene
rowData(ods)\$counted1sample = rowSums(assay(ods)) > 0

# External data check
if (is.null(ods@colData\$GENE_COUNTS_FILE)){ #column does not exist in sample annotation table
    has_external <- FALSE
}else if(all(is.na(ods@colData\$GENE_COUNTS_FILE))){ #column exists but it has no values
    has_external <- FALSE
}else if(all(ods@colData\$GENE_COUNTS_FILE == "")){ #column exists with non-NA values but this group has all empty strings
    has_external <- FALSE
}else{ #column exists with non-NA values and this group has at least 1 non-empty string
    has_external <- TRUE
}

if(has_external){
    ods@colData\$isExternal <- as.factor(!(ods@colData\$GENE_COUNTS_FILE %in% c("", NA)))
}else{
    ods@colData\$isExternal <- as.factor(FALSE)
}


# Save the ods before filtering to preserve the original number of genes
saveRDS(ods, paste(prefix, ".Rds", sep=""))


## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    bioconductor-genomicfeatures:', as.character(packageVersion('GenomicFeatures'))),
        paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
        paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
    ),
'versions.yml')
