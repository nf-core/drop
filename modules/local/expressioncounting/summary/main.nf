process EXPRESSIONCOUNTING_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c76d0f6d35c84c2e9f9b883f1d7683c1fbdc09e8b8019f38b0915cc038dfd76f/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-outrider_bioconductor-summarizedexperiment_r-base_pruned:beedc80a74b3909a' }"

    input:
    tuple val(meta), path(ods)
    val(bam_cov)

    output:
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Summary.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- file("${prefix}.Rds", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-cowplot:', as.character(packageVersion('cowplot'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
            paste('    r-ggthemes:', as.character(packageVersion('ggthemes'))),
            paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
            paste('    bioconductor-GenomicAlignments::', as.character(packageVersion('GenomicAlignments'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
