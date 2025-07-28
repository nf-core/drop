#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_4_countRNA_nonSplitReads_merge.R
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_5_countRNA_collect.R
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/03_filter_expression_FraseR.R
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/03_filter_expression_FraseR.R
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/04_fit_hyperparameters_FraseR.R
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/FRASER/05_fit_autoencoder_FraseR.R

### COUNTRNA_NONSPLITREADSMERGE

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

dataset    <- "$drop_group"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Move non split counts to the correct cache directory
system("mv cache/nonSplicedCounts/raw-local cache/nonSplicedCounts/raw-local-${drop_group}")

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS("$non_split_counts_granges")

# Remove the cache file of merged counts from previous runs
nonSplitCountsDir <- file.path(workingDir, "savedObjects",
                            paste0("raw-local-", dataset), 'nonSplitCounts')
if(dir.exists(nonSplitCountsDir)){
    unlink(nonSplitCountsDir, recursive = TRUE)
}

# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds,
                                                    splitCountRanges=splitCounts_gRanges,
                                                    minAnchor=5,
                                                    recount=FALSE,
                                                    longRead=${long_read ? "TRUE" : "FALSE"},)

# remove cache dir as its not needed later
if(dir.exists(nonSplitCountsDir)){
    unlink(nonSplitCountsDir, recursive = TRUE)
}

message(date(), ":", dataset, " nonSplit counts done")

### COUNTRNA_COLLECT

dataset    <- "$drop_group"
workingDir <- "./"
saveDir    <- "savedObjects/raw-local-$drop_group"

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))
splitCounts_gRanges <- readRDS("$splitcounts_granges")
spliceSiteCoords <- readRDS("$splicesites_splitcounts")

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts_h5 <- HDF5Array::HDF5Array(file.path(saveDir, "rawCountsJ.h5"), "rawCountsJ")
splitCounts_se <- SummarizedExperiment(
    colData = colData(fds),
    rowRanges = splitCounts_gRanges,
    assays = list(rawCountsJ=splitCounts_h5)
)


nonSplitCounts_h5 <- HDF5Array::HDF5Array(file.path(saveDir, "rawCountsSS.h5"), "rawCountsSS")
nonSplitCounts_se <- SummarizedExperiment(
    colData = colData(fds),
    rowRanges = spliceSiteCoords,
    assays = list(rawCountsSS=nonSplitCounts_h5)
)

# Add Counts to FRASER object
fds <- addCountsToFraserDataSet(fds=fds, splitCounts=splitCounts_se,
                                nonSplitCounts=nonSplitCounts_se)

# Save final FRASER object
fds <- saveFraserDataSet(fds)

### FRASER_PSIVALUECALCULATION

dataset    <- "$drop_group"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Calculating PSI values
fds <- calculatePSIValues(fds)

# FRASER object after PSI value calculation
fds <- saveFraserDataSet(fds)

### FRASER_FILTEREXPRESSION

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

### FRASER_FITHYPERPARAMETERS

rseed <- ${random_seed instanceof Integer ? random_seed : random_seed ? 'TRUE' : 'FALSE'}
if (!isFALSE(rseed)) {
    if(isTRUE(rseed)){
        set.seed(42)
    } else if (is.numeric(rseed)){
        set.seed(as.integer(rseed))
    }
}

#+ input
dataset    <- "$drop_group"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Load PSI data
fds <- loadFraserDataSet(dir=workingDir, name=dataset)
fitMetrics(fds) <- psiTypes

# Run hyper parameter optimization
implementation <- "$implementation"
mp <- $max_tested_dimension_proportion

# Get range for latent space dimension
a <- 2
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

maxSteps <- 12
if(mp < 6){
    maxSteps <- 15
}

Nsteps <- min(maxSteps, b)
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique

for(type in psiTypes){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type,
                            implementation=implementation,
                            q_param=pars_q,
                            plot = FALSE)
    fds <- saveFraserDataSet(fds)
}
fds <- saveFraserDataSet(fds)

### FRASER_FITAUTOENCODER

dataset    <- "$drop_group"
workingDir <- "./"

register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

fds <- loadFraserDataSet(dir=workingDir, name=dataset)

# Fit autoencoder
# run it for every type
implementation <- "$implementation"

for(type in psiTypes){
    currentType(fds) <- type
    q <- bestQ(fds, type)
    verbose(fds) <- 3   # Add verbosity to the FRASER object
    fds <- fit(fds, q=q, type=type, iterations=15, implementation=implementation)
    fds <- saveFraserDataSet(fds)
}

# remove .h5 files from previous runs with other FRASER version
fdsDir <- "$fds/$drop_group/"
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
