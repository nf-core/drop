process FILTERCOUNTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a971e624aad46fc672638427d28931d13257642ce54743d6a3869d305eecd642/data' :
        'community.wave.seqera.io/library/bioconductor-genomicfeatures_bioconductor-outrider_bioconductor-summarizedexperiment_r-base_r-data.table:c6ad593656300741' }"

    input:
    tuple val(meta), path(counts), path(txdb)
    val(fpkmCutoff)

    output:
    tuple val(meta), path("*.Rds")  , emit: output
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'filterCounts.R'

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
            paste('    bioconductor-genomicfeatures:', as.character(packageVersion('GenomicFeatures'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
