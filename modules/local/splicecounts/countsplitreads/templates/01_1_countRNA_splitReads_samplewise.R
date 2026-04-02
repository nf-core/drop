#!/usr/bin/env Rscript
# https://github.com/gagneurlab/drop/blob/master/drop/modules/aberrant-splicing-pipeline/Counting/01_1_countRNA_splitReads_samplewise.R

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
library(BSgenome)

dataset    <- "$drop_group"
genomeAssembly <- "$genome_assembly"

# Read FRASER object
fds <- loadFraserDataSet(dir="./", name=paste0("raw-local-", dataset))

# Get sample id from wildcard
sample_id <- "$sample_id"

# If data is not strand specific, add genome info
genome <- NULL
if(strandSpecific(fds[, sample_id]) == 0){
    genome <- getBSgenome(genomeAssembly)
}

# Count splitReads for a given sample id
sample_result <- countSplitReads(sampleID = sample_id,
                                fds = fds,
                                NcpuPerSample = ${task.cpus},
                                recount = ${recount ? 'TRUE' : 'FALSE'},
                                keepNonStandardChromosomes = ${keep_non_standard_chrs ? 'TRUE' : 'FALSE'},
                                genome = genome)

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (split counts) = ", length(sample_result))

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
