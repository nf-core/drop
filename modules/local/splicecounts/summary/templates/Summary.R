#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/Summary.R

source("$config", echo=FALSE)
configure_fraser("$fraser_version")

suppressPackageStartupMessages({
    library(cowplot)
})

opts_chunk\$set(fig.width=10, fig.height=8)

#+ input
dataset    <- "$drop_group"
workingDir <- "./"

fdsLocal <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))
fdsMerge <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

has_external <- !(all(is.na(fdsMerge@colData\$SPLICE_COUNTS_DIR)) || is.null(fdsMerge@colData\$SPLICE_COUNTS_DIR))
if(has_external){
    fdsMerge@colData\$isExternal <- as.factor(!is.na(fdsMerge@colData\$SPLICE_COUNTS_DIR))
}else{
    fdsMerge@colData\$isExternal <- as.factor(FALSE)
}
devNull <- saveFraserDataSet(fdsMerge,dir=workingDir, name=paste0("raw-", dataset))


#' ## Number of samples:
#' Local: `r sum(!as.logical(fdsMerge@colData\$isExternal))`
#' External: `r sum(as.logical(fdsMerge@colData\$isExternal))`
n_local <- sum(!as.logical(fdsMerge@colData\$isExternal))
n_external <- sum(as.logical(fdsMerge@colData\$isExternal))
count_df <- data.frame(sample_type = c("Local", "External"),
    count = c(n_local, n_external), stringsAsFactors = FALSE)
write.table(count_df, file="number_of_samples_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#' ## Number of introns:
#' Local (before filtering): `r length(rowRanges(fdsLocal, type = "j"))`
#' ```{asis, echo = has_external}
#' Merged with external counts (before filtering):
#' ```
#' ```{r, eval = has_external, echo=FALSE}
#' length(rowRanges(fdsMerge, type = "j"))
#' ```
#' After filtering: `r sum(mcols(fdsMerge, type="j")[,"passed"])`
n_introns_local_before  <- length(rowRanges(fdsLocal, type = "j"))
n_introns_merged_before <- if (has_external) length(rowRanges(fdsMerge, type = "j")) else NA_integer_
n_introns_after         <- sum(mcols(fdsMerge, type="j")[,"passed"])
introns_df <- data.frame(type = c("Local (before filtering)", if (has_external) "Merged with external counts (before filtering)" else NULL, "After filtering"),
    count = c(n_introns_local_before, if (has_external) n_introns_merged_before else NULL, n_introns_after), stringsAsFactors = FALSE)
write.table(introns_df, file="number_of_introns_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#' ## Number of splice sites:
#' Local: `r length(rowRanges(fdsLocal, type = "theta"))`
#' ```{asis, echo = has_external}
#' Merged with external counts:
#' ```
#' ```{r, eval = has_external, echo=FALSE}
#' length(rowRanges(fdsMerge, type = "theta"))
#' ```
n_splice_local  <- length(rowRanges(fdsLocal, type = "theta"))
n_splice_merged <- if (has_external) length(rowRanges(fdsMerge, type = "theta")) else NA_integer_
splice_df <- data.frame(type = c("Local", if (has_external) "Merged with external counts" else NULL),
    count = c(n_splice_local, if (has_external) n_splice_merged else NULL), stringsAsFactors = FALSE)
write.table(splice_df, file="number_of_splice_sites_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


#' ```{asis, echo = has_external}
#' ## Comparison of local and external counts
#' **Using external counts**
#' External counts introduce some complexity into the problem of counting junctions
#' because it is unknown whether or not a junction is not counted (because there are no reads)
#' compared to filtered and not present due to legal/personal sharing reasons. As a result,
#' after merging the local (counted from BAM files) counts and the external counts, only the junctions that are
#' present in both remain. As a result it is likely that the number of junctions will decrease after merging.
#'
#' ```
#' ```{r, eval = has_external, echo=has_external}
if(has_external){
    externalCountIDs <- colData(fdsMerge)[as.logical(colData(fdsMerge)[,"isExternal"]),"sampleID"]
    localCountIDs <- colData(fdsMerge)[!as.logical(colData(fdsMerge)[,"isExternal"]),"sampleID"]
    cts <- K(fdsMerge,"psi5")
    ctsLocal<- cts[,localCountIDs,drop=FALSE]
    ctsExt<- cts[,externalCountIDs,drop=FALSE]
    rowMeanLocal <- rowMeans(ctsLocal)
    rowMeanExt <- rowMeans(ctsExt)
    dt <- data.table("Mean counts of local samples" = rowMeanLocal,
                     "Mean counts of external samples" = rowMeanExt)
    p <- ggplot(dt,aes(x = `Mean counts of local samples`, y= `Mean counts of external samples`)) +
       geom_hex() + theme_cowplot(font_size = 16) +
	   theme_bw() + scale_x_log10() + scale_y_log10() +
       geom_abline(slope = 1, intercept =0) +
       scale_color_brewer(palette="Dark2")
    ggsave("comparison_local_and_external_mqc.png", p, width = 5, height = 4, dpi = 196, bg = "white")
}
#' ```
#'

#' ## Expression filtering
#' The expression filtering step removes introns that are lowly expressed. The requirements for an intron to pass this filter are:
#'
#' * at least 1 sample has `r snakemake@config\$aberrantSplicing\$minExpressionInOneSample` counts (K) for the intron
#' * at least `r 100*(1-snakemake@config\$aberrantSplicing\$quantileForFiltering)`% of the samples need to have a total of at least `r snakemake@config\$aberrantSplicing\$quantileMinExpression` reads for the splice metric denominator (N) of the intron
png("expression_filtering_mqc.png", width = 8, height = 6.4, units = "in", res = 196)
plotFilterExpression(fdsMerge) +
    labs(title="", x="Mean Intron Expression", y="Introns") +
    theme_cowplot(font_size = 16)
dev.off()

#' ## Variability filtering
#' The variability filtering step removes introns that have no or little variability in the splice metric values across samples. The requirement for an intron to pass this filter is:
#'
#' * at least 1 sample has a difference of at least `r snakemake@config\$aberrantSplicing\$minDeltaPsi` in the splice metric compared to the mean splice metric of the intron
png("variability_filtering_mqc.png", width = 8, height = 6.4, units = "in", res = 196)
plotFilterVariability(fdsMerge) +
    labs(title="", y="Introns") +
    theme_cowplot(font_size = 16)
dev.off()

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
        paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
        paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
        paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
        paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
    ),
'versions.yml')
