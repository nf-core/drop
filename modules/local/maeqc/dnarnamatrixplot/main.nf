process MAEQC_DNARNAMATRIXPLOT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8f75e4dbdb0d66ca82fb85cf1e41c02cc80f7f4f70d53af77c37b8bc815b7f3e/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-ggplot2_r-pheatmap_pruned:fd75147f7ce60ef1' }"

    input:
    tuple val(meta) , path(mat_qc), path(samplesheet), val(drop_group)
    val(dnaRnaMatchCutoff)

    output:
    tuple val(meta), path("matching_values_distribution_mqc.png")  , emit: matching_values_distribution
    tuple val(meta), path("heatmap_matching_variants_mqc.*")       , emit: heatmap_matching_variants
    tuple val(meta), path("identify_matching_samples_mqc.tsv")     , emit: identify_matching_samples
    tuple val(meta), path("false_matches_mqc.tsv")                 , emit: false_matches
    tuple val(meta), path("false_mismatches_mqc.tsv")              , emit: false_mismatches
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'DNA_RNA_matrix_plot.R'

    stub:
    """
    #!/usr/bin/env Rscript
    a <- file("matching_values_distribution_mqc.png", "w")
    close(a)
    a <- file("heatmap_matching_variants_mqc.png", "w")
    close(a)
    a <- file("identify_matching_samples_mqc.tsv", "w")
    close(a)
    a <- file("false_matches_mqc.tsv", "w")
    close(a)
    a <- file("false_mismatches_mqc.tsv", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-reshape2:', as.character(packageVersion('reshape2'))),
            paste('    r-pheatmap:', as.character(packageVersion('pheatmap'))),
            paste('    r-rcolorbrewer:', as.character(packageVersion('RColorBrewer')))
        ),
    'versions.yml')
    """
}
