process SPLICECOUNTS_MERGENONSPLITREADS {
    tag "${meta.id}"
    label 'process_low'

    // TODO discuss if this container should be split up, library loading should be done in the module script instead in that case
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ffa8a0cabe0d92e58c583162eb9fa946edcffc5a16958405f070c99ce93a7c59/data' :
        'community.wave.seqera.io/library/bioconductor-bsgenome.hsapiens.ncbi.grch38_bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-bsgenome_pruned:0776613eaf49b161' }"

    input:
    tuple val(meta), path(fds, stageAs: "input/savedObjects/*"), path(non_split_counts_granges), path(non_split_reads, stageAs:"input/cache/nonSplicedCounts/raw-local/*"), path(bams), path(bais), val(drop_group)
    val(long_read)
    val(recount)
    val(fraser_version)

    output:
    tuple val(meta), path("cache/nonSplicedCounts/raw-local-${drop_group}"), emit: cache
    tuple val(meta), path("savedObjects/raw-local-${drop_group}" )         , emit: fdsobj
    path  "versions.yml"                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template '01_4_countRNA_nonSplitReads_merge.R'

    stub:
    """
    #!/usr/bin/env Rscript
    dir.create("cache/nonSplicedCounts/raw-local-${drop_group}", recursive = TRUE)
    a <- file("cache/nonSplicedCounts/raw-local-${drop_group}/spliceSitesCoordinates.RDS", "w")
    close(a)
    dir.create("savedObjects/raw-local-${drop_group}", recursive = TRUE)
    a <- file("savedObjects/raw-local-${drop_group}/rawCountsSS.h5", "w")
    close(a)
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
    """
}
