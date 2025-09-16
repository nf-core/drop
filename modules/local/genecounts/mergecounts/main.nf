process GENECOUNTS_MERGECOUNTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/eca2d19aed7918a344093e4627135c21403e08bb16b730d924cd93ae5b428464/data' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-summarizedexperiment_r-base_r-data.table_pruned:c96d6e13d480e41d' }"

    input:
    tuple val(meta), path(counts), path(count_ranges), val(external_counts_ids)
    tuple val(meta2), path(samplesheet) // Revisit this later, should be a samplesheet in CSV or TSV format

    output:
    tuple val(meta), path("*.Rds")   , emit: output
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'mergeCounts.R'

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
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-data.table:', as.character(packageVersion('dplyr'))),
            paste('    bioconductor-biocparallel:', as.character(packageVersion('BiocParallel'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment')))
        ),
    'versions.yml')
    """
}
