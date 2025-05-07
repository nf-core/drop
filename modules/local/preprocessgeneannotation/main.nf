process PREPROCESSGENEANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3c/3c54f46cc4a1f9f92c216095525003cce5ad6581696086a76d9df45767d45a48/data' :
        'community.wave.seqera.io/library/bioconductor-genomicfeatures_bioconductor-rtracklayer_bioconductor-txdbmaker_r-base_pruned:2d189ec410925335' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.db")                       , emit: txdb
    tuple val(meta), path("*.count_ranges.Rds")         , emit: count_ranges
    tuple val(meta), path("*.gene_name_mapping.tsv")    , emit: gene_name_mapping
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'preprocessGeneAnnotation.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(GenomicFeatures)
        library(data.table)
        library(rtracklayer)
        library(magrittr)
    })
    a <- file("${prefix}.db", "w")
    close(a)
    a <- file("${prefix}.count_ranges.Rds", "w")
    close(a)
    a <- file("${prefix}.gene_name_mapping.tsv", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    bioconductor-genomicfeatures:', as.character(packageVersion('GenomicFeatures'))),
            paste('    bioconductor-rtracklayer:', as.character(packageVersion('rtracklayer'))),
            paste('    bioconductor-txdbmaker:', as.character(packageVersion('txdbmaker')))
        ),
    'versions.yml')
    """
}
