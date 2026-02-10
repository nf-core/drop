process SPLICECOUNTS_MERGESPLITREADS {
    tag "${meta.id}"
    label 'process_low'

    // TODO discuss if this container should be split up, library loading should be done in the module script instead in that case
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/4319818c9416d02ea4e4bd81bad32533ce5268b81a5c53b7117984edad6655f1/data' :
        'community.wave.seqera.io/library/bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-bsgenome_bioconductor-delayedmatrixstats_pruned:e8af96bc80316fdc' }"

    input:
    tuple val(meta), path(fds, stageAs: "input/savedObjects/*"), path(split_counts, stageAs: "input/cache/splitCounts/*"), path(bams), path(bais), val(drop_group)
    val(min_expression_in_one_sample)
    val(recount)
    val(fraser_version)
    path(config) // Pass "${projectDir}/assets/helpers/aberrant_splicing_config.R" to this input

    output:
    tuple val(meta), path("cache/raw-local-${drop_group}")       , emit: cache
    tuple val(meta), path("savedObjects/raw-local-${drop_group}"), emit: fdsobj
    path  "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template '01_2_countRNA_splitReads_merge.R'

    stub:
    """
    #!/usr/bin/env Rscript
    dir.create("cache/raw-local-${drop_group}", recursive = TRUE)
    a <- file("cache/raw-local-${drop_group}/gRanges_splitCounts.rds", "w")
    close(a)
    a <- file("cache/raw-local-${drop_group}/gRanges_NonSplitCounts.rds", "w")
    close(a)
    a <- file("cache/raw-local-${drop_group}/spliceSites_splitCounts.rds", "w")
    close(a)
    dir.create("savedObjects/raw-local-${drop_group}", recursive = TRUE)
    a <- file("savedObjects/raw-local-${drop_group}/rawCountsJ.h5", "w")
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
