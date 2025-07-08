#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/mae-pipeline/QC/create_matrix_dna_rna_cor.R

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(GenomicRanges)
    library(magrittr)
    library(BiocParallel)
    library(data.table)
})

register(MulticoreParam($task.cpus))
sa <- fread("$samplesheet",
            colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))

# Read the test vcf as GRanges
gr_test <- readVcf("$qc_vcf") %>% granges()
mcols(gr_test)\$GT <- "0/0"

# Obtain the rna and vcf files
rna_samples <- as.character(${rna_ids.collect { "\"$it\"" }.join(', ')})
mae_res <- c(${res.collect { "\"$it\"" }.join(', ')})

rows_in_group <- sapply(strsplit(sa\$DROP_GROUP, ',|, '), function(d) "$drop_group" %in% d)
vcf_cols <- sa[rows_in_group, .(DNA_ID, DNA_VCF_FILE)] %>% unique
dna_samples <- vcf_cols\$DNA_ID
vcf_files <- c(${vcfs.collect { "\"$it\"" }.join(', ')})

# Read all RNA genotypes into a list
rna_gt_list <- bplapply(mae_res, readRDS)

# For every DNA, find the overlapping variants with each RNA
N <- length(vcf_files)

lp <- bplapply(1:N, function(i){

    # Read sample vcf file

    sample <- dna_samples[i] %>% as.character()

    param <-  ScanVcfParam(fixed=NA, info='NT', geno='GT', samples=sample, trimEmpty=TRUE)
    vcf_sample <- readVcf(vcf_files[i], param = param, row.names = FALSE)
    # Get GRanges and add Genotype
    gr_sample <- granges(vcf_sample)

    if(!is.null(geno(vcf_sample)\$GT)){
        gt <- geno(vcf_sample)\$GT
        gt <- gsub('0|0', '0/0', gt, fixed = TRUE)
        gt <- gsub('0|1', '0/1', gt, fixed = TRUE)
        gt <- gsub('1|0', '0/1', gt, fixed = TRUE)
        gt <- gsub('1|1', '1/1', gt, fixed = TRUE)
    } else if(!is.null(info(vcf_sample)\$NT)){
        gt <- info(vcf_sample)\$NT
        gt <- gsub('ref', '0/0', gt)
        gt <- gsub('het', '0/1', gt)
        gt <- gsub('hom', '1/1', gt)
    }

    mcols(gr_sample)\$GT <- gt

    # Find overlaps between test and sample
    gr_res <- copy(gr_test)
    seqlevelsStyle(gr_res) <- seqlevelsStyle(seqlevelsInUse(gr_sample))[1]
    ov <- findOverlaps(gr_res, gr_sample, type = 'equal')
    mcols(gr_res)[from(ov),]\$GT <- mcols(gr_sample)[to(ov),]\$GT

    # Find similarity between DNA sample and RNA sample
    x <- vapply(rna_gt_list, function(gr_rna){
        seqlevelsStyle(gr_rna) <- seqlevelsStyle(gr_res)[1]
        ov <- findOverlaps(gr_res, gr_rna, type = 'equal')
        gt_dna <- gr_res[from(ov)]\$GT
        gt_rna <- gr_rna[to(ov)]\$RNA_GT
        sum(gt_dna == gt_rna) / length(gt_dna)
    }, 1.0)
    return(x)
})

# Create a matrix
mat <- do.call(rbind, lp)
row.names(mat) <- dna_samples
colnames(mat) <- rna_samples
# mat <- mat[sa[rows_in_group, DNA_ID], sa[rows_in_group, RNA_ID],drop=FALSE]

saveRDS(mat, 'dna_rna_qc_matrix.Rds')

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    bioconductor-biocparallel:', as.character(packageVersion('BiocParallel'))),
        paste('    bioconductor-variantannotation:', as.character(packageVersion('VariantAnnotation'))),
        paste('    bioconductor-genomicranges:', as.character(packageVersion('GenomicRanges')))
    ),
'versions.yml')
