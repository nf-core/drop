process OUTRIDER_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02b44404c5bf450fc2ec0a3247408b09096e0f429afc4c4ad4cf46f06238a110/data' :
        'community.wave.seqera.io/library/bioconductor-outrider_bioconductor-summarizedexperiment_r-base_r-cowplot_pruned:0059eff556e3bc37' }"

    input:
    tuple val(meta), path(ods), path(results), val(drop_group), val(annotation_id)
    val(padjCutoff)
    val(zScoreCutoff)

    output:
    tuple val(meta), path("outrider_overview_mqc.tsv")             , emit: sample_mqc
    tuple val(meta), path("Encoding_dimension_mqc.png")            , emit: encoding_mqc
    tuple val(meta), path("Aberrant_genes_per_sample_mqc.png")     , emit: genes_mqc
    tuple val(meta), path("Batch_correction_raw_mqc.png")          , emit: batch_mqc_1
    tuple val(meta), path("Batch_correction_normalized_mqc.png")   , emit: batch_mqc_2
    tuple val(meta), path("geneSampleHeatmap_raw_mqc.png")         , emit: heatmap_mqc_1
    tuple val(meta), path("geneSampleHeatmap_normalized_mqc.png")  , emit: heatmap_mqc_2
    tuple val(meta), path("BCV_mqc.png")                           , emit: bcv_mqc
    tuple val(meta), path("outrider_result_overview_mqc.tsv")      , emit: result_mqc_1
    tuple val(meta), path("Aberrant_samples_mqc.tsv")              , emit: result_mqc_2
    tuple val(meta), path("significant_results_mqc.tsv")           , emit: result_mqc_3
    tuple val(meta), path("OUTRIDER_results_*.tsv")                , emit: result
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Summary.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- file("outrider_overview_mqc.tsv", "w")
    close(a)
    a <- file("Encoding_dimension_mqc.png", "w")
    close(a)
    a <- file("Aberrant_genes_per_sample_mqc.png", "w")
    close(a)
    a <- file("Batch_correction_raw_mqc.png", "w")
    close(a)
    a <- file("Batch_correction_normalized_mqc.png", "w")
    close(a)
    a <- file("geneSampleHeatmap_raw_mqc.png", "w")
    close(a)
    a <- file("geneSampleHeatmap_normalized_mqc.png", "w")
    close(a)
    a <- file("BCV_mqc.png", "w")
    close(a)
    a <- file("outrider_result_overview_mqc.tsv", "w")
    close(a)
    a <- file("Aberrant_samples_mqc.tsv", "w")
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
