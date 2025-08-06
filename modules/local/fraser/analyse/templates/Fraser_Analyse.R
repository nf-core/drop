#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/06_annotate_genes.R
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/07_calculation_stats_FraseR.R
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/08_extract_results_FraseR.R

### FRASER_ANNOTATEGENES

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

annotation <- "$annotation_id"
dataset    <- "$drop_group"
fdsFile    <- "$fds"
workingDir <- "./"
outputDir  <- "./processed_data/"

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

### FRASER_CALCULATIONSTATS

annotation <- "$annotation_id"
dataset    <- "$drop_group"
fdsFile    <- "$fds"
workingDir <- "./processed_data/"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# read in subsets from sample anno if present (returns NULL if not present)
source("$parse_subsets_for_FDR")
fraser_sample_ids <- list(${sample_ids.collect { "\"$it\"" }.join(', ')})
subsets <- parse_subsets_for_FDR("$genes_to_test",
                                sampleIDs=fraser_sample_ids)

# Load FRASER data
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# Calculate stats
for (type in psiTypes) {
    # Pvalues
    fds <- calculatePvalues(fds, type=type)
    # Adjust Pvalues
    fds <- calculatePadjValues(fds, type=type, subsets=subsets)
}

fds <- saveFraserDataSet(fds, dir="./")

# remove .h5 files from previous runs with other FRASER version
fdsDir <- "$fds/${drop_group}--${annotation_id}/"
pvalFiles <- grep("p(.*)BetaBinomial_(.*).h5",
                list.files(fdsDir),
                value=TRUE)
for(type in psiTypesNotUsed){
    pvalFilesType <- grep(type, pvalFiles, value=TRUE)
    for(pFile in pvalFilesType){
        unlink(file.path(fdsDir, pFile))
    }
}

### FRASER_EXTRACTRESULTS

source("$add_HPO_cols")
library(AnnotationDbi)

annotation    <- "$annotation_id"
dataset    <- "$drop_group"
fdsFile    <- "$fds"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Load fds and create a new one
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# Extract results per junction
res_junc <- results(fds, psiType=psiTypes,
                    padjCutoff=$padj_cutoff,
                    deltaPsiCutoff=$delta_psi_cutoff,)
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')

# Add features
if(nrow(res_junc_dt) > 0){
    # number of samples per gene and variant
    res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
    res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
} else{
    warning("The aberrant splicing pipeline gave 0 intron-level results for the ", dataset, " dataset.")
}

# Extract full results by gene
res_gene <- results(fds, psiType=psiTypes,
                    aggregate=TRUE, collapse=FALSE,
                    all=TRUE)
res_genes_dt   <- as.data.table(res_gene)
print('Results per gene extracted')
write_tsv(res_genes_dt, file="results_gene_all.tsv")

# Subset gene results to aberrant
padj_cols <- grep("padjustGene", colnames(res_genes_dt), value=TRUE)
res_genes_dt <- res_genes_dt[do.call(pmin, c(res_genes_dt[,padj_cols, with=FALSE],
                                                list(na.rm = TRUE))) <= $padj_cutoff &
                                    abs(deltaPsi) >= $delta_psi_cutoff &
                                    totalCounts >= 5,]

if(nrow(res_genes_dt) > 0){
    # add HPO overlap information
    sa <- fread("$samplesheet",
                colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
    if(!is.null(sa\$HPO_TERMS)){
        if(!all(is.na(sa\$HPO_TERMS)) & ! all(sa\$HPO_TERMS == '')){
            res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = "$hpo_file")
        }
    }
} else{
    warning("The aberrant splicing pipeline gave 0 gene-level results for the ", dataset, " dataset.")
}

# Annotate results with spliceEventType and blacklist region overlap
txdb <- loadDb("$txdb")

# annotate the type of splice event and UTR overlap
if(nrow(res_junc_dt) > 0){
    res_junc_dt <- annotatePotentialImpact(result=res_junc_dt, txdb=txdb, fds=fds)
}
if(nrow(res_genes_dt) > 0){
    res_genes_dt <- annotatePotentialImpact(result=res_genes_dt, txdb=txdb, fds=fds)
}

# set genome assembly version to load correct blacklist region BED file (hg19 or hg38)
assemblyVersion <- "$genome_assembly"
if(grepl("grch37", assemblyVersion, ignore.case=TRUE)) assemblyVersion <- "hg19"
if(grepl("grch38", assemblyVersion, ignore.case=TRUE)) assemblyVersion <- "hg38"

# annotate overlap with blacklist regions
if(assemblyVersion %in% c("hg19", "hg38")){
    if(nrow(res_junc_dt) > 0){
        res_junc_dt <- flagBlacklistRegions(result=res_junc_dt,
                                        assemblyVersion=assemblyVersion)
    }
    if(nrow(res_genes_dt) > 0){
        res_genes_dt <- flagBlacklistRegions(result=res_genes_dt,
                                        assemblyVersion=assemblyVersion)
    }
} else{
    message(date(), ": cannot annotate blacklist regions as no blacklist region\n",
            "BED file is available for genome assembly version ", assemblyVersion,
            " as part of FRASER.")
}

# Results
write_tsv(res_junc_dt, file="results_per_junction.tsv")
write_tsv(res_genes_dt, file="results.tsv")

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
