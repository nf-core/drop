#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-expression-pipeline/Counting/countReads.R

suppressPackageStartupMessages({
    library(data.table)
    library(Rsamtools)
    library(BiocParallel)
    library(GenomicAlignments)
})

# Initialize the prefix from the module
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

# Rename the BAM file to make sure the correct ID is used
if ("$bam" != "${meta.id}.bam") {
    file.rename("$bam", "${meta.id}.bam")
}

# Get strand specific information from sample annotation

strand <- tolower("$strand")
count_mode <- "$count_mode"
paired_end <- as.logical("$paired_end")
overlap <- as.logical("$count_overlaps")
inter_feature <- ! overlap # inter_feature = FALSE does not allow overlaps

# infer preprocessing and strand info
preprocess_reads <- NULL
if (strand == "yes") {
    strand_spec <- T
} else if (strand == "no"){
    strand_spec <- F
} else if (strand == "reverse") {
    # set preprocess function for later
    preprocess_reads <- invertStrand
    strand_spec <- T
} else {
    stop(paste("invalid strand information", strand))
}

# read files
bam_file <- BamFile("${meta.id}.bam", yieldSize = $yield_size)
count_ranges <- readRDS("$count_ranges")
# set chromosome style
seqlevelsStyle(count_ranges) <- seqlevelsStyle(bam_file)

# show info
message(paste("input:", "${meta.id}.bam"))
message(paste("output:", paste(prefix, ".Rds", sep="")))
message(paste('\tcount mode:', count_mode, sep = "\t"))
message(paste('\tpaired end:', paired_end, sep = "\t"))
message(paste('\tinter.feature:', inter_feature, sep = "\t"))
message(paste('\tstrand:', strand, sep = "\t"))
message(paste('\tstrand specific:', strand_spec, sep = "\t"))
message(paste(seqlevels(count_ranges), collapse = ' '))

# start counting
message("\ncounting")
starttime <- Sys.time()
se <- summarizeOverlaps(
    count_ranges
    , bam_file
    , mode = count_mode
    , singleEnd = !paired_end
    , ignore.strand = !strand_spec  # FALSE if done strand specifically
    , fragments = F
    , count.mapped.reads = T
    , inter.feature = inter_feature # TRUE: reads mapping to multiple features are dropped
    , preprocess.reads = preprocess_reads
    , BPPARAM = MulticoreParam($task.cpus)
)
saveRDS(se, paste(prefix, ".Rds", sep=""))
message("done")

print(format(Sys.time()- starttime)) # print time taken
print(sum(assay(se))) # total counts

# Write the versions file
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
        paste('    bioconductor-rsamtools:', as.character(packageVersion('Rsamtools'))),
        paste('    bioconductor-biocparallel:', as.character(packageVersion('BiocParallel')))
    ),
'versions.yml')
