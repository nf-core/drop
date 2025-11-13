process OUTRIDER_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/22/2234ffd3a93c24242604e355c216fda3154e3522d8d80fadf7c11451c7229d41/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_bioconductor-outrider_bioconductor-summarizedexperiment_pruned:a73c412fd56d7e1c' }"

    input:
    tuple val(meta), path(ods), path(results), val(drop_group), val(annotation_id)
    val(padjCutoff)
    val(zScoreCutoff)

    output:
    tuple val(meta), path("outrider_overview_mqc.tsv")             , emit: outrider_overview
    tuple val(meta), path("encoding_dimension_mqc.png")            , emit: encoding_dimension
    tuple val(meta), path("aberrant_genes_per_sample_mqc.png")     , emit: aberrant_genes_per_sample
    tuple val(meta), path("batch_correction_raw_mqc.png")          , emit: batch_correction_raw
    tuple val(meta), path("batch_correction_normalized_mqc.png")   , emit: batch_correction_normalized
    tuple val(meta), path("geneSampleHeatmap_raw_mqc.png")         , emit: geneSampleHeatmap_raw
    tuple val(meta), path("geneSampleHeatmap_normalized_mqc.png")  , emit: geneSampleHeatmap_normalized
    tuple val(meta), path("bcv_mqc.png")                           , emit: bcv
    tuple val(meta), path("outrider_result_overview_mqc.tsv")      , emit: outrider_result_overview
    tuple val(meta), path("aberrant_samples_mqc.tsv")              , emit: aberrant_samples
    tuple val(meta), path("significant_results_mqc.tsv")           , emit: significant_results
    tuple val(meta), path("OUTRIDER_results_*.tsv")                , emit: results
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Summary.R'

    stub:
    """
    #!/usr/bin/env Rscript
    a <- file("outrider_overview_mqc.tsv", "w")
    close(a)
    a <- file("encoding_dimension_mqc.png", "w")
    close(a)
    a <- file("aberrant_genes_per_sample_mqc.png", "w")
    close(a)
    a <- file("batch_correction_raw_mqc.png", "w")
    close(a)
    a <- file("batch_correction_normalized_mqc.png", "w")
    close(a)
    a <- file("geneSampleHeatmap_raw_mqc.png", "w")
    close(a)
    a <- file("geneSampleHeatmap_normalized_mqc.png", "w")
    close(a)
    a <- file("bcv_mqc.png", "w")
    close(a)
    a <- file("outrider_result_overview_mqc.tsv", "w")
    close(a)
    a <- file("aberrant_samples_mqc.tsv", "w")
    close(a)
    a <- file("significant_results_mqc.tsv", "w")
    close(a)
    a <- file("OUTRIDER_results_.tsv", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-cowplot:', as.character(packageVersion('cowplot'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
            paste('    r-ggthemes:', as.character(packageVersion('ggthemes'))),
            paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
