#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/exportCounts.R

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

library(AnnotationDbi)

#
# input
#
annotation_file <- "${annotation}"
workingDir <- "."
dataset    <- "${drop_group}"

out_k_files <- c("k_j_counts.tsv.gz", "k_theta_counts.tsv.gz")
out_n_files <- c("n_psi5_counts.tsv.gz", "n_psi3_counts.tsv.gz", "n_theta_counts.tsv.gz")

# Read annotation and extract known junctions
txdb <- loadDb(annotation_file)
introns <- unique(unlist(intronsByTranscript(txdb)))
introns <- keepStandardChromosomes(introns, pruning.mode = 'coarse')
length(introns)

# Read FRASER object, adapt chr style and subset to known junctions
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))
seqlevels(fds) <- seqlevelsInUse(fds)
seqlevelsStyle(fds) <- seqlevelsStyle(introns)[1]

# recalculate theta.h5 in Jaccard metric
if (fitMetrics(fds) == "jaccard") {
    fds <- calculatePSIValues(fds, type="theta")
} else {
    fds <- calculatePSIValues(fds, type="jaccard")
}

fds_known <- fds[unique(to(findOverlaps(introns, rowRanges(fds, type="j"), type="equal"))),]

# save k/n counts
sapply(c(out_k_files, out_n_files), function(i){
    ctsType <- toupper(strsplit(basename(i), "_")[[1]][1])
    psiType <- strsplit(basename(i), "_")[[1]][2]

    cts <- as.data.table(get(ctsType)(fds_known, type=psiType))
    grAnno <- rowRanges(fds_known, type=psiType)
    anno <- as.data.table(grAnno)
    anno <- anno[,.(seqnames, start, end, strand)]

    fwrite(cbind(anno, cts), file=i, quote=FALSE, row.names=FALSE, sep="\t", compress="gzip")
}) %>% invisible()

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
        paste('    bioconductor-annotationdbi:', as.character(packageVersion('AnnotationDbi'))),
        paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
        paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
        paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
        paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
    ),
'versions.yml')
