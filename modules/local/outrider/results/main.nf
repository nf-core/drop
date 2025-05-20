process OUTRIDER_RESULTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c0ec1df146f59c0d2fa6b535330c2a0b133f1533f9ec6caca6a2a445b16e5791/data' :
        'community.wave.seqera.io/library/bioconductor-outrider_bioconductor-summarizedexperiment_r-base_r-data.table_pruned:242dc589aa3ec32c' }"

    input:
    tuple val(meta), path(ods)
    tuple val(meta2), path(gene_name_mapping)
    tuple val(meta3), path(samplesheet) // revisit this later, should be a file in CSV or TSV format
    tuple val(meta4), path(hpo_file)
    val(padjCutoff)
    val(zScoreCutoff)
    path(add_HPO_cols) // Pass ${projectDir}/assets/helpers/add_HPO_cols.R to this input

    output:
    tuple val(meta), path("*.tsv")  , emit: results
    tuple val(meta), path("*.Rds")  , emit: results_all
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'results.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    a <- file("${prefix}.Rds", "w")
    close(a)
    a <- file("${prefix}.tsv", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
            paste('    bioconductor-summarizedexperiments:', as.character(packageVersion('SummarizedExperiment'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
