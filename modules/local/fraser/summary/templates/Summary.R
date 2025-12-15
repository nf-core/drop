#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/Summary.R

source("$config", echo=FALSE)
if("$fraser_version" == "FRASER2"){
    pseudocount(0.1)
    psiTypes <- c("jaccard")
    psiTypesNotUsed <- c("psi5", "psi3", "theta")
} else{
    pseudocount(1)
    psiTypes <- c("psi5", "psi3", "theta")
    psiTypesNotUsed <- c("jaccard")
}

suppressPackageStartupMessages({
  library(cowplot)
  library(RColorBrewer)
})

#+ input
dataset    <- "$drop_group"
annotation <- "$annotation_id"
padj_cutoff <- $padj_cutoff
deltaPsi_cutoff <- $deltaPsi_cutoff

fds <- loadFraserDataSet(dir="./", name=paste(dataset, annotation, sep = '--'))
hasExternal <- length(levels(colData(fds)\$isExternal) > 1)

#' Number of samples: `r nrow(colData(fds))`
#'
#' Number of introns: `r length(rowRanges(fds, type = "j"))`
#'
#' Number of splice sites: `r length(rowRanges(fds, type = "ss"))`
n_samples <- nrow(colData(fds))
n_introns <- length(rowRanges(fds, type = "j"))
n_splice <- length(rowRanges(fds, type = "ss"))

count_df <- data.frame(sample_type = c("Number of samples", "Number of introns","Number of splice sites"),
    count = c(n_samples, n_introns, n_splice), stringsAsFactors = FALSE)
write.table(count_df, file="fraser_overview_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# used for most plots
dataset_title <- paste0("Dataset: ", dataset, "--", annotation)


#' ## Estimation of the optimal latent space dimension
if (isTRUE(metadata(fds)[["useOHTtoObtainQ"]])){
  type <- "jaccard"
  g <- plotEncDimSearch(fds, type=type, plotType="sv")
  g <- g + theme_cowplot(font_size = 16) +
    ggtitle(paste0("Q estimation, ", type)) + theme(legend.position = "none")
  png(paste0("q_estimation_", type,"_mqc.png"), width = 4, height = 3, units = "in", res = 196)
  print(g)
  dev.off()
} else{
  for(type in psiTypes){
    g <- plotEncDimSearch(fds, type=type, plotType="auc")
    g <- g + theme_cowplot(font_size = 16) +
      ggtitle(paste0("Q estimation, ", type)) + theme(legend.position = "none")
    png(paste0("q_estimation_", type,"_mqc.png"), width = 4, height = 3, units = "in", res = 196)
    print(g)
    dev.off()
  }
}

#' ## Aberrantly spliced genes per sample
png("aberrantly_spliced_genes_mqc.png", width = 4, height = 3, units = "in", res = 196)
plotAberrantPerSample(fds, type=psiTypes,
                      padjCutoff = padj_cutoff, deltaPsiCutoff = deltaPsi_cutoff,
                      aggregate=TRUE, main=dataset_title) +
  theme_cowplot(font_size = 16) +
  theme(legend.position = "top")
dev.off()

#' ## Batch Correlation: samples x samples
topN <- 30000
topJ <- 10000
anno_color_scheme <- brewer.pal(n = 3, name = 'Dark2')[1:2]

for(type in psiTypes){
  for(normalized in c(F,T)){
    hm <- plotCountCorHeatmap(
      object=fds,
      type = type,
      logit = TRUE,
      topN = topN,
      topJ = topJ,
      plotType = "sampleCorrelation",
      normalized = normalized,
      annotation_col = "isExternal",
      annotation_row = NA,
      sampleCluster = NA,
      minDeltaPsi = minDeltaPsi,
      plotMeanPsi = FALSE,
      plotCov = FALSE,
      annotation_legend = TRUE,
      annotation_colors = list(isExternal = c("FALSE" = anno_color_scheme[1], "TRUE" = anno_color_scheme[2]))
    )
    ggsave(
      filename = sprintf("batch_correlation_%s_%s_mqc.png", type, normalized),
      plot     = hm,
      width    = 5,
      height   = 4,
      units    = "in",
      dpi      = 196
    )
    dev.off()
  }
}

#' ## Results
res <- fread("$results")

#' Total gene-level splicing outliers: `r nrow(res)`
#'
n_outliers <- nrow(res)
count_df <- data.frame(type = c("Total gene-level splicing outliers"),
    count = c(n_outliers), stringsAsFactors = FALSE)
write.table(count_df, file="total_outliers_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# file <- snakemake@output\$res_html
# write_tsv(res, file = file)
# #+ echo=FALSE, results='asis'
# cat(paste0("<a href='./", basename(file), "'>Download FRASER results table</a>"))

# round numbers
if(nrow(res) > 0){
  res[, pValue := signif(pValue, 3)]
  res[, deltaPsi := signif(deltaPsi, 2)]
  res[, psiValue := signif(psiValue, 2)]
  res[, pValueGene := signif(pValueGene, 2)]
  padjGene_cols <- grep("padjustGene", colnames(res), value=TRUE)
  for(padj_col in padjGene_cols){
      res[, c(padj_col) := signif(get(padj_col), 2)]
  }
}

# DT::datatable(
#   head(res, 1000),
#   caption = 'FRASER results (up to 1,000 rows shown)',
#   options=list(scrollX=TRUE),
#   escape=FALSE,
#   filter = 'top'
# )
setorder(res, pValue)
fwrite(head(cbind(row_id = rownames(res), res), 500),"results_mqc.tsv", sep = "\t", quote = F)

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-rmarkdown:', as.character(packageVersion('rmarkdown'))),
        paste('    r-knitr:', as.character(packageVersion('knitr'))),
        paste('    r-devtools:', as.character(packageVersion('devtools'))),
        paste('    r-yaml:', as.character(packageVersion('yaml'))),
        paste('    r-bbmisc:', as.character(packageVersion('BBmisc'))),
        paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
        paste('    r-plotly:', as.character(packageVersion('plotly'))),
        paste('    r-rhdf5:', as.character(packageVersion('rhdf5'))),
        paste('    r-cowplot:', as.character(packageVersion('cowplot'))),
        paste('    r-rcolorbrewer:', as.character(packageVersion('RColorBrewer'))),
        paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
        paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
        paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
        paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
    ),
'versions.yml')
