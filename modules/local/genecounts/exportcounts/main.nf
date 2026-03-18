process GENECOUNTS_EXPORTCOUNTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed29816b0f1418e3c91cfaeb3ef53006ec05f012766ddd13956799fc057ef8f2/data' :
        'community.wave.seqera.io/library/bioconductor-outrider_bioconductor-summarizedexperiment_r-base_r-data.table:7bac555dbe5a5fc5' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.tsv.gz")   , emit: counts
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'exportCounts.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- gzfile("${prefix}.tsv.gz", "w")
    close(a)

    # Write the versions file
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
