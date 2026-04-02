#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/03_filter_expression_FraseR.R

source(Sys.which("aberrant_splicing_config.R"), echo=FALSE)
configure_fraser("$fraser_version")

opts_chunk\$set(fig.width=12, fig.height=8)

# input
dataset    <- "$drop_group"
workingDir <- "./"
exCountIDs <- c(${external_count_ids.collect { id -> "\"$id\"" }.join(', ')})
exCountFiles <- c(${splice_counts_dirs.collect { file -> "\"$file\"" }.join(', ')})
sample_anno_file <- "$samplesheet"
minExpressionInOneSample <- $min_expressions_in_one_sample
quantile <- $quantile_filtering
quantileMinExpression <- $quantile_min_expression
minDeltaPsi <- $min_delta_psi
filterOnJaccard <- ("$fraser_version" == "FRASER2")

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Add external data if provided by dataset
if(length(exCountIDs) > 0){
    message("create new merged fraser object")
    fds <- saveFraserDataSet(fds,dir = workingDir, name=paste0("raw-", dataset))

    for(resource_raw in unique(exCountFiles)){
        # NEXTFLOW: untar external splice counts if necessary
        resource <- resource_raw
        if(endsWith(resource_raw, 'tar.gz')) {
            untar(resource_raw)
            resource <- file.path(workingDir, sub("\\\\.tar\\\\.gz\$", "", resource_raw))
        }
        exSampleIDs <- exCountIDs[exCountFiles == resource_raw]
        exAnno <- fread(sample_anno_file, key="RNA_ID")[J(exSampleIDs)]
        exAnno[, strand:=exAnno\$STRAND]
        exAnno\$strand<- sapply(exAnno[, strand], switch, 'no' = 0L, 'unstranded' = 0L,
                                    'yes' = 1L, 'stranded' = 1L, 'reverse' = 2L, -1L)
        setnames(exAnno, "RNA_ID", "sampleID")

        ctsNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
        ctsFiles <- paste0(resource, "/", ctsNames, "_counts.tsv.gz")

        # Merging external counts restricts the junctions to those that
        # are only present in both the counted (fromBam) junctions AND the
        # junctions from the external counts.
        fds <- mergeExternalData(fds=fds, countFiles=ctsFiles,
                sampleIDs=exSampleIDs, annotation=exAnno)
        fds@colData\$isExternal <- as.factor(!is.na(fds@colData\$SPLICE_COUNTS_DIR))
    }
} else {
    fds@colData\$isExternal <- as.factor(FALSE)
    workingDir(fds) <- workingDir
    name(fds) <- paste0("raw-", dataset)
}

# filter for expression and write it out to disc.
fds <- filterExpressionAndVariability(fds,
        minExpressionInOneSample = minExpressionInOneSample,
        quantile=quantile,
        quantileMinExpression=quantileMinExpression,
        minDeltaPsi = minDeltaPsi,
        filterOnJaccard=filterOnJaccard,
        filter=FALSE)

devNull <- saveFraserDataSet(fds,dir = workingDir)

# Keep junctions that pass filter
name(fds) <- dataset
if (${filter ? 'TRUE' : 'FALSE'} == TRUE) {
    filtered <- mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}

seqlevels(fds) <- seqlevelsInUse(fds)
colData(fds)\$sampleID <- as.character(colData(fds)\$sampleID)
fds <- saveFraserDataSet(fds,dir = workingDir)

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
