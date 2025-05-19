#!/usr/bin/env Rscript
#https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-expression-pipeline/OUTRIDER/results.R

source("$add_HPO_cols")

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

ods <- readRDS("$ods")
res <- results(ods, padjCutoff = $padjCutoff,
            zScoreCutoff = $zScoreCutoff, all = TRUE)

# Add fold change
res[, foldChange := round(2^l2fc, 2)]

# Save all the results and significant ones
saveRDS(res, paste(prefix, ".Rds", sep=""))

# Subset to significant results
padj_cols <- grep("padjust", colnames(res), value=TRUE)
res <- res[do.call(pmin, c(res[,padj_cols, with=FALSE], list(na.rm = TRUE)))
                <= $padjCutoff &
            abs(zScore) >= $zScoreCutoff]

gene_annot_dt <- fread("$gene_name_mapping")
if(!is.null(gene_annot_dt\$gene_name)){
    if(grepl('ENSG00', res[1,geneID]) & grepl('ENSG00', gene_annot_dt[1,gene_id])){
        res <- merge(res, gene_annot_dt[, .(gene_id, gene_name)],
                    by.x = 'geneID', by.y = 'gene_id', sort = FALSE, all.x = TRUE)
        setnames(res, 'gene_name', 'hgncSymbol')
        res <- cbind(res[, .(hgncSymbol)], res[, - 'hgncSymbol'])
    }
}

# Add HPO terms, requires online connection and for there to be annotated HPO terms
# TODO rethink this in a more Nextflow way
sa <- fread("$samplesheet",
                colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
if(!is.null(sa\$HPO_TERMS) & nrow(res) > 0){
    if(!all(is.na(sa\$HPO_TERMS)) & ! all(sa\$HPO_TERMS == '')){
        res <- add_HPO_cols(res, hpo_file = "$hpo_file")
    }
}


# Save results
fwrite(res, paste(prefix, ".tsv", sep=""), sep = "\t", quote = F)

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
        paste('    bioconductor-summarizedexperiments:', as.character(packageVersion('SummarizedExperiment'))),
        paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
    ),
'versions.yml')
