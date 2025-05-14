process MERGECOUNTS {
    tag "${meta.id}"
    label 'process_low'

    // TODO: Add singularity container once seqera containers works again
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-summarizedexperiment_r-base_r-data.table_r-dplyr:a612bfc9c2b3630e' }"

    input:
    // TODO: Filter out the exclude ID files in the subworkflow
    tuple val(meta), path(counts)
    tuple val(meta2), path(count_ranges)
    tuple val(meta3), path(samplesheet) // Revisit this later, should be a samplesheet in CSV or TSV format
    val(exclude_ids)

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
