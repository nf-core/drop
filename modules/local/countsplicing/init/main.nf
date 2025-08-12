process COUNTSPLICING_INIT {
    tag "${meta.id}"
    label 'process_single'

    // TODO discuss if this container should be split up, library loading should be done in the module script instead in that case
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc5b31e03713c86445cb08c02014dda9298c2c0a7b20a17cda43abc66f6eb63c/data' :
        'community.wave.seqera.io/library/bioconductor-delayedmatrixstats_bioconductor-fraser_bioconductor-genomicalignments_r-base_pruned:8f3a9483f1399b0a' }"

    input:
    tuple val(meta), path(col_data), val(drop_group)
    val(fraser_version)
    path(config) // Pass "${projectDir}/assets/helpers/aberrant_splicing_config.R" to this input

    output:
    tuple val(meta), path("savedObjects/raw-local-${drop_group}") , emit: fdsobj
    path  "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template '01_0_countRNA_init.R'

    stub:
    """
    #!/usr/bin/env Rscript
    dir.create("savedObjects")
    dir.create("savedObjects/raw-local-${drop_group}")
    a <- file("savedObjects/raw-local-${drop_group}/fds-object.RDS", "w")
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
            paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
        ),
    'versions.yml')
    """
}
