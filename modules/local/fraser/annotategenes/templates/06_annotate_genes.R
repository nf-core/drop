#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/06_annotate_genes.R

source(Sys.which("aberrant_splicing_config.R"), echo=FALSE)
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

annotation <- "$annotation_id"
dataset    <- "$drop_group"
fdsFile    <- "$fds"
workingDir <- "./"
outputDir  <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Load fds and create a new one
fds_input <- loadFraserDataSet(dir=workingDir, name=dataset)

# Read annotation and match the chr style
txdb <- loadDb("$txdb")
orgdb <- fread("$gene_name_mapping")

seqlevels_fds <- seqlevelsStyle(fds_input)[1]
seqlevelsStyle(orgdb\$seqnames) <- seqlevels_fds
seqlevelsStyle(txdb) <- seqlevels_fds

# Annotate the fds with gene names and save it as a new object
fds_input <- annotateRangesWithTxDb(fds_input, txdb = txdb, orgDb = orgdb,
                    feature = 'gene_name', featureName = 'hgnc_symbol',
                    keytype = 'gene_id')

# add basic annotations for overlap with the reference annotation
# run this function before creating the results table to include it there
fds_input <- annotateIntronReferenceOverlap(fds_input, txdb)

# save fds
fds <- saveFraserDataSet(fds_input, dir=outputDir, name = paste(dataset, annotation, sep = '--'), rewrite = TRUE)

# remove .h5 files from previous runs with other FRASER version
fdsDir <- "savedObjects/${drop_group}--${annotation_id}"
for(type in psiTypesNotUsed){
    predMeansFile <- file.path(fdsDir, paste0("predictedMeans_", type, ".h5"))
    if(file.exists(predMeansFile)){
        unlink(predMeansFile)
    }
}

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
        paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
        paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
        paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
        paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
    ),
'versions.yml')
