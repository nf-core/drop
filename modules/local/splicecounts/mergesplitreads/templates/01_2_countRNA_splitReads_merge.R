#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_2_countRNA_splitReads_merge.R

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
workingDir <- "."
minExpressionInOneSample <- $min_expression_in_one_sample


register(MulticoreParam($task.cpus))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam($task.cpus))

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# If samples are recounted, remove the merged ones
splitCountsDir <- file.path(workingDir, "savedObjects",
                            paste0("raw-local-", dataset), 'splitCounts')
if(${recount ? 'TRUE': 'FALSE'} == TRUE & dir.exists(splitCountsDir)){
    unlink(splitCountsDir, recursive = TRUE)
}

# Get and merge splitReads for all sample ids
splitCounts <- getSplitReadCountsForAllSamples(fds=fds,
                                            recount=FALSE)
# Extract, annotate and save granges
splitCountRanges <- rowRanges(splitCounts)

dir.create(paste0("cache/raw-local-", dataset))

# Annotate granges from the split counts
splitCountRanges <- FRASER:::annotateSpliceSite(splitCountRanges)
saveRDS(splitCountRanges, paste0("cache/raw-local-", dataset, "/gRanges_splitCounts.rds"))
# additionally save as tsv.gz (for easier AbSplice input)
fwrite(as.data.table(splitCountRanges),
        gsub(".Rds", ".tsv.gz", paste0("cache/raw-local-", dataset, "/gRanges_splitCounts.rds"),
            ignore.case=TRUE))

# Create ranges for non split counts
# Subset by minExpression
maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= minExpressionInOneSample
# extract granges after filtering
splitCountRanges <- splitCountRanges[passed,]

saveRDS(splitCountRanges, paste0("cache/raw-local-", dataset, "/gRanges_NonSplitCounts.rds"))

# Extract splitSiteCoodinates: extract donor and acceptor sites
# take either filtered or full fds
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges)
saveRDS(spliceSiteCoords, paste0("cache/raw-local-", dataset, "/spliceSites_splitCounts.rds"))


message(date(), ": ", dataset, " total no. splice junctions = ",
        length(splitCounts))

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
