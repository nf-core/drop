#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_2_countRNA_splitReads_merge.R

source("$config", echo=FALSE)

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

# Annotate granges from the split counts
splitCountRanges <- FRASER:::annotateSpliceSite(splitCountRanges)
saveRDS(splitCountRanges, paste("cache/raw-local-", dataset, "/gRanges_splitCounts.rds"))
# additionally save as tsv.gz (for easier AbSplice input)
fwrite(as.data.table(splitCountRanges),
        gsub(".Rds", ".tsv.gz", paste("cache/raw-local-", dataset, "/gRanges_splitCounts.rds"),
            ignore.case=TRUE))

# Create ranges for non split counts
# Subset by minExpression
maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= minExpressionInOneSample
# extract granges after filtering
splitCountRanges <- splitCountRanges[passed,]

saveRDS(splitCountRanges, paste("cache/raw-local-", dataset, "/gRanges_NonSplitCounts.rds"))

# Extract splitSiteCoodinates: extract donor and acceptor sites
# take either filtered or full fds
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges)
saveRDS(spliceSiteCoords, paste("cache/raw-local-", dataset, "/spliceSites_splitCounts.rds"))


message(date(), ": ", dataset, " total no. splice junctions = ",
        length(splitCounts))
