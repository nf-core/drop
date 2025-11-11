process OUTRIDER_FIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_bioconductor-outrider_bioconductor-summarizedexperiment_pruned:f837c8409d1d5f5c' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_bioconductor-outrider_bioconductor-summarizedexperiment_pruned:a73c412fd56d7e1c' }"

    input:
    tuple val(meta), path(ods)
    val(implementation)
    val(ae_max_tested_dimension_proportion)
    val(random_seed)
    val(ae_use_grid_search_to_obtain_q)

    output:
    tuple val(meta), path("*.Rds")            , emit: ods_fitted
    path  "versions.yml"                      , emit: versions

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
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
