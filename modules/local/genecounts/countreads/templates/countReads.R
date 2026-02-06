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

# Rename the BAM file and its index to make sure the correct ID is used
 if ("$bam" != "${meta.id}.bam") {
    file.rename("$bam", "${meta.id}.bam")
    # Also rename the BAM index file
    bam_index <- paste0("$bam", ".bai")
    if (file.exists(bam_index)) {
        file.rename(bam_index, "${meta.id}.bam.bai")
        # Update timestamp to be newer than BAM file
        Sys.setFileTime("${meta.id}.bam.bai", Sys.time())
    }
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
bam_file <- BamFile("$bam", yieldSize = $yield_size)
count_ranges <- readRDS("$count_ranges")
# set chromosome style
seqlevelsStyle(count_ranges) <- seqlevelsStyle(bam_file)

# run it in parallel across all chromosomes
gene_seqnames = as.character(sapply(seqnames(count_ranges), runValue))

# show info
message(paste("input:", "$bam"))
message(paste("output:", paste(prefix, ".Rds", sep="")))
message(paste('\tcount mode:', count_mode, sep = "\t"))
message(paste('\tpaired end:', paired_end, sep = "\t"))
message(paste('\tinter.feature:', inter_feature, sep = "\t"))
message(paste('\tstrand:', strand, sep = "\t"))
message(paste('\tstrand specific:', strand_spec, sep = "\t"))
message(paste(seqlevels(count_ranges), collapse = ' '))

# function to count per chromosome
count_per_chromosome <- function(i, gene_seqnames, count_ranges, bam_file, yield_size,
                                 count_mode, paired_end, strand_spec, preprocess_reads) {
  genes_to_count <- gene_seqnames == i
  message(paste("counting seqname '", i, "' with '", sum(genes_to_count), "' genes.", sep = ""))

  # subset the count_ranges for the current seqname
  current_ranges <- count_ranges[genes_to_count]

  # create a new BamFile for the current seqname
  current_bam_file <- BamFile("$bam", yieldSize=$yield_size)
  # summarize overlaps
  summarizeOverlaps(
    current_ranges,
    current_bam_file,
    mode = count_mode,
    singleEnd = !paired_end,
    ignore.strand = !strand_spec,
    fragments = FALSE,
    count.mapped.reads = TRUE,
    preprocess.reads = preprocess_reads,
    inter.feature = inter_feature,
    param = ScanBamParam(which = GRanges(i, IRanges(1, seqlengths(bam_file)[i])))
  )
}

# run counting in parallel across all chromosomes
message("\nstarting counting expression ...")

starttime <- Sys.time()

iterate <- seqlevels(count_ranges)
bpparam <- SerialParam()
counts <- bplapply(iterate, count_per_chromosome,
          gene_seqnames, count_ranges, bam_file, $yield_size,
          count_mode, paired_end, strand_spec, preprocess_reads,
          BPPARAM = bpparam
)
# merge SE objects - concatenate by rows (genes) across chromosomes
se <- do.call(rbind, counts)
colnames(se) <- "{$meta.id}"

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
