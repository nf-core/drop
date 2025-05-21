process OUTRIDER_RUN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/6586ceea612bc211a55c10d444e660c29be0ca3087538a06fce79267d07a82f7/data' :
        'community.wave.seqera.io/library/bioconductor-outrider_bioconductor-summarizedexperiment_r-base_r-dplyr_pruned:0b2146d44ef9304b' }"

    input:
    tuple val(meta), path(ods)

    output:
    tuple val(meta), path("*fitted.Rds")            , emit: ods_fitted
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'runOutrider.R'

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
            paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
            paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
            paste('    r-tools:', as.character(packageVersion('tools'))),
            paste('    r-yaml:', as.character(packageVersion('yaml'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
