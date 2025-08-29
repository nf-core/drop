#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/mae-pipeline/MAE/Results.R

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(tidyr)
    library(GenomicRanges)
    library(SummarizedExperiment)
    library(R.utils)
    library(dplyr)
})

# Read all MAE results files
rmae <- lapply(c(${mae_res.collect { "\"$it\"" }.join(', ')}), function(m){
    rt <- readRDS(m)
    # force consistant UCSC chromosome style
    rt <- rt[!grepl("chr",contig),contig:= paste0("chr",contig)]
    return(rt)
}) %>% rbindlist()

# re-factor contig
rmae\$contig <- factor(rmae\$contig)

# Convert results into GRanges
rmae_ranges <- GRanges(seqnames = rmae\$contig,
                        IRanges(start = rmae\$position, end = rmae\$position), strand = '*')

# Read annotation and convert into GRanges
gene_annot_dt <- fread("$gene_name_mapping")
gene_annot_ranges <- GRanges(seqnames = gene_annot_dt\$seqnames,
                            IRanges(start = gene_annot_dt\$start, end = gene_annot_dt\$end),
                            strand = gene_annot_dt\$strand)
gene_annot_ranges <- keepStandardChromosomes(gene_annot_ranges, pruning.mode = 'coarse')

# Keep the chr style of the annotation in case the results contain different styles
seqlevelsStyle(rmae_ranges) <- seqlevelsStyle(gene_annot_ranges)

# Overlap results and annotation
fo <- findOverlaps(rmae_ranges, gene_annot_ranges)

# Add the gene names
res_annot <- cbind(rmae[from(fo), ],  gene_annot_dt[to(fo), .(gene_name, gene_type)])

# Prioritize protein coding genes
res_annot <- rbind(res_annot[gene_type == 'protein_coding'],
                    res_annot[gene_type != 'protein_coding'])

# Write all the other genes in another column
res_annot[, aux := paste(contig, position, sep = "-")]
rvar <- unique(res_annot[, .(aux, gene_name)])
rvar[, N := 1:.N, by = aux]

r_other <- rvar[N > 1, .(other_names = paste(gene_name, collapse = ',')), by = aux]
res <- merge(res_annot, r_other, by = 'aux', sort = FALSE, all.x = TRUE)
res[, c('aux') := NULL]
res <- res[, .SD[1], by = .(ID, contig, position)]

# Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

# Calculate variant frequency within cohort
maxCohortFreq <- $max_cohort_freq
res[, N_var := .N, by = .(gene_name, contig, position)]
res[, cohort_freq := round(N_var / uniqueN(ID), 3)]

res[, rare := (rare | is.na(rare)) & cohort_freq <= maxCohortFreq]

# Add significance columns
allelicRatioCutoff <- $allelic_ratio_cutoff
res[, MAE := padj <= $padj_cutoff &
        (altRatio >= allelicRatioCutoff | altRatio <= (1-allelicRatioCutoff))
    ]
res[, MAE_ALT := MAE == TRUE & altRatio >= allelicRatioCutoff]

n_samples <- uniqueN(res\$ID)
n_genes <- uniqueN(res\$gene_name)
n_significant <- uniqueN(res[MAE_ALT == TRUE, ID])
count_df <- data.frame(sample_type = c("Number of samples", "Number of genes", "Number of samples with significant MAE for ALT"),
    count = c(n_samples, n_genes, n_significant), stringsAsFactors = FALSE)
write.table(count_df, file="mae_overview_mqc.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Save full results zipped
res[, altRatio := round(altRatio, 3)]
fwrite(res, paste0("MAE_results_all_", prefix, ".tsv.gz"), sep = '\t',
        row.names = F, quote = F, compress = 'gzip')

# Save significant results
fwrite(res[MAE_ALT == TRUE], paste0("MAE_results_", prefix, ".tsv"), sep = '\t', row.names = F, quote = F)

# Save significant results
fwrite(res[MAE_ALT == TRUE & rare == TRUE], paste0("MAE_results_", prefix, "_rare.tsv"), sep = '\t', row.names = F, quote = F)


# Add columns for plot
res[, N := .N, by = ID]
plot_res <- res[,.(N = .N,
            N_MAE = sum(MAE==T),
            N_MAE_REF=sum(MAE==T & MAE_ALT == F),
            N_MAE_ALT=sum(MAE_ALT == T),
            N_MAE_REF_RARE = sum(MAE ==T & MAE_ALT==F & rare == T),
            N_MAE_ALT_RARE = sum(MAE_ALT ==T & rare ==T)
            ),by = ID]


melt_dt <- melt(plot_res, id.vars = 'ID')
melt_dt[variable == 'N', variable := '>10 counts']
melt_dt[variable == 'N_MAE', variable := '+MAE']
melt_dt[variable == 'N_MAE_REF', variable := '+MAE for\nREF']
melt_dt[variable == 'N_MAE_ALT', variable := '+MAE for\nALT']
melt_dt[variable == 'N_MAE_REF_RARE', variable := '+MAE for REF\n& rare']
melt_dt[variable == 'N_MAE_ALT_RARE', variable := '+MAE for ALT\n& rare']

## Cascade plot
p <- ggplot(melt_dt, aes(variable, value)) + geom_boxplot() +
    scale_y_log10(limits = c(1,NA)) + theme_bw(base_size = 14) +
    labs(y = 'Heterozygous SNVs per patient', x = '') +
        annotation_logticks(sides = "l")
ggsave("cascade_plot_mqc.png", p, width = 5, height = 4, dpi = 196, bg = "white")

#'
#' ## Variant Frequency within Cohort
p <- ggplot(unique(res[,cohort_freq,by =.(gene_name, contig, position)]),aes(x = cohort_freq)) + geom_histogram( binwidth = 0.02)  +
    geom_vline(xintercept = maxCohortFreq, col = "red",linetype="dashed") + theme_bw(base_size = 14) +
    xlim(0,NA) + xlab("Variant frequency in cohort") + ylab("Variants")
ggsave("variant_frequency_mqc.png", p, width = 5, height = 4, dpi = 196, bg = "white")

#' Median of each category
# DT::datatable(melt_dt[, .(median = median(value, na.rm = T)), by = variable])
fwrite(melt_dt[, .(median = median(value, na.rm = T)), by = variable],"median_of_each_category_mqc.tsv", sep = "\t", quote = F)

# round numbers
if(nrow(res) > 0){
    res[, pvalue := signif(pvalue, 3)]
    res[, padj := signif(padj, 3)]
    res[, log2FC := signif(log2FC, 3)]
}
#'
#' ## MAE Results table
# DT::datatable(
#     head(res[MAE_ALT == TRUE], 1000),
#     caption = 'MAE results (up to 1,000 rows shown)',
#     options=list(scrollX=TRUE),
#     filter = 'top'
# )
res_alt <- res[MAE_ALT == TRUE]
fwrite(head(cbind(row_id = rownames(res_alt), res_alt), 500),"results_mqc.tsv", sep = "\t", quote = F)

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    r-r.utils:', as.character(packageVersion('R.utils'))),
        paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
        paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
        paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
        paste('    r-dt:', as.character(packageVersion('DT'))),
        paste('    bioconductor-genomicranges:', as.character(packageVersion('GenomicRanges'))),
        paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment')))
    ),
'versions.yml')
