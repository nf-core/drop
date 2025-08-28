#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-expression-pipeline/OUTRIDER/Summary.R

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(cowplot)
    library(data.table)
    library(dplyr)
    library(ggthemes)
})

# used for most plots
dataset_title <- paste("Dataset:", paste("$drop_group", "$annotation_id", sep = "--"))

ods <- readRDS("$ods")
if (is.null(colData(ods)\$isExternal)) colData(ods)\$isExternal <- FALSE

#' Number of samples: `r ncol(ods)`
#' Number of expressed genes: `r nrow(ods)`
n_samples <- ncol(ods)
n_expressed_genes <- nrow(ods)
count_df <- data.frame(sample_type = c("Number of samples", "Number of expressed genes"),
    count = c(n_samples, n_expressed_genes), stringsAsFactors = FALSE)
write.table(count_df, file="outrider_overview_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#'
#' ## Visualize
#' ### Encoding dimension
p <- plotEncDimSearch(ods) + labs(title = dataset_title) + theme_cowplot() + background_grid() +
    scale_color_brewer(palette = "Dark2")
ggsave("encoding_dimension_mqc.png", p, width = 5, height = 3.75, dpi = 196, bg = "white")

#' ### Aberrantly expressed genes per sample
p <- plotAberrantPerSample(ods, main = dataset_title, padjCutoff = $padjCutoff,
    zScoreCutoff = $zScoreCutoff)
ggsave("aberrant_genes_per_sample_mqc.png", p, width = 5, height = 3.75, dpi = 196,
    bg = "white")

#' ### Batch correction
#+ countCorHeatmap, fig.height=8, fig.width=8
png("batch_correction_raw_mqc.png", width = 5, height = 5, units = "in", res = 196)
plotCountCorHeatmap(ods, normalized = FALSE, colGroups = "isExternal", colColSet = "Dark2",
    main = paste0("Raw Counts (", dataset_title, ")"))
dev.off()
png("batch_correction_normalized_mqc.png", width = 5, height = 5, units = "in", res = 196)
plotCountCorHeatmap(ods, normalized = TRUE, colGroups = "isExternal", colColSet = "Dark2",
    main = paste0("Normalized Counts (", dataset_title, ")"))
dev.off()

#' ### Expression by gene per sample
#+ geneSampleHeatmap, fig.height=12, fig.width=8
png("geneSampleHeatmap_raw_mqc.png", width = 6, height = 9, units = "in", res = 196)
plotCountGeneSampleHeatmap(ods, normalized = FALSE, nGenes = 50, colGroups = "isExternal",
    colColSet = "Dark2", main = paste0("Raw Counts (", dataset_title, ")"), bcvQuantile = 0.95,
    show_names = "row")
dev.off()
png("geneSampleHeatmap_normalized_mqc.png", width = 6, height = 9, units = "in",
    res = 196)
plotCountGeneSampleHeatmap(ods, normalized = TRUE, nGenes = 50, colGroups = "isExternal",
    colColSet = "Dark2", main = paste0("Normalized Counts (", dataset_title, ")"),
    bcvQuantile = 0.95, show_names = "row")
dev.off()

#' ### BCV - Biological coefficient of variation
# function to calculate BCV before autoencoder
estimateThetaWithoutAutoCorrect <- function(ods) {

    ods1 <- OutriderDataSet(countData = counts(ods), colData = colData(ods))
    # use rowMeans as expected means
    normalizationFactors(ods1) <- matrix(rowMeans(counts(ods1)), ncol = ncol(ods1),
        nrow = nrow(ods1))
    ods1 <- fit(ods1)
    theta(ods1)

    return(theta(ods1))
}

before <- data.table(when = "Before", BCV = 1/sqrt(estimateThetaWithoutAutoCorrect(ods)))
after <- data.table(when = "After", BCV = 1/sqrt(theta(ods)))
bcv_dt <- rbind(before, after)

# boxplot of BCV Before and After Autoencoder
#+ BCV, fig.height=5, fig.width=6
p <- ggplot(bcv_dt, aes(when, BCV)) + geom_boxplot() + theme_bw(base_size = 14) +
    labs(x = "Autoencoder correction", y = "Biological coefficient \nof variation",
        title = dataset_title)
ggsave("bcv_mqc.png", p, width = 3, height = 3, dpi = 196, bg = "white")

#' ## Results
res <- fread("$results")

#' Total number of expression outliers: `r nrow(res)`
#' Samples with at least one outlier gene: `r res[, uniqueN(sampleID)]`
n_outliers <- nrow(res)
n_samples_outliers <- res[, uniqueN(sampleID)]
count_df <- data.frame(sample_type = c("Total number of expression outliers", "Samples with at least one outlier gene"),
    count = c(n_outliers, n_samples_outliers), stringsAsFactors = FALSE)
write.table(count_df, file="outrider_result_overview_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#'
#' ### Aberrant samples
#'
#' An aberrant sample is one that has more than 0.1% of the total genes called as outliers.
if (nrow(res) > 0) {
    ab_table <- res[AberrantBySample > nrow(ods)/1000, .(`Outlier genes` = .N), by = .(sampleID)] %>%
        unique
    if (nrow(ab_table) > 0) {
        setorder(ab_table, "Outlier genes")
        writeLines(c('# plot_type: "table"', '# format: "tsv"'), con = file.path("aberrant_samples_mqc.tsv"))
        write.table(ab_table, file="aberrant_samples_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE, sep = "\t")
    } else {
        writeLines(c('# plot_type: "html"','no aberrant samples.'), con = file.path("aberrant_samples_mqc.tsv"))
    }
} else {
    writeLines(c('# plot_type: "html"','no aberrant samples.'), con = file.path("aberrant_samples_mqc.tsv"))
}


#' ## Results table

## Save results table in the html folder and provide link to download
file <- paste0("OUTRIDER_results_", "$drop_group", ".tsv")
fwrite(res, file, sep = "\t", quote = F)


if (nrow(res) > 0) {
    res[, pValue := format(pValue, scientific = T, digits = 3)]
    res[, padjust := format(padjust, scientific = T, digits = 3)]

    # DT::datatable(head(res, 1000), caption = 'OUTRIDER results (up to 1,000
    # rows shown)', options=list(scrollX=TRUE), filter = 'top')
    writeLines(c('# plot_type: "table"', '# format: "tsv"'), con = file.path("significant_results_mqc.tsv"))
    write.table(head(cbind(row_id = rownames(res), res), 100), file="significant_results_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE, sep = "\t")
} else {
    writeLines(c('# plot_type: "html"','no significant results.'), con = file.path("significant_results_mqc.tsv"))
}

## VERSIONS FILE
writeLines(c("\"${task.process}\":", paste("    r-base:", strsplit(version[["version.string"]],
    " ")[[1]][3]), paste("    r-cowplot:", as.character(packageVersion("cowplot"))),
    paste("    r-data.table:", as.character(packageVersion("data.table"))), paste("    r-ggplot2:",
        as.character(packageVersion("ggplot2"))), paste("    r-ggthemes:", as.character(packageVersion("ggthemes"))),
    paste("    r-dplyr:", as.character(packageVersion("dplyr"))), paste("    bioconductor-summarizedexperiment:",
        as.character(packageVersion("SummarizedExperiment"))), paste("    bioconductor-outrider:",
        as.character(packageVersion("OUTRIDER")))), "versions.yml")
