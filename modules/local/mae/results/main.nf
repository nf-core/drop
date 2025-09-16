process MAE_RESULTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e08488ec876636112342a762061843003bfbf25f8e66ae2925951d2cc8304992/data' :
        'community.wave.seqera.io/library/bioconductor-genomicranges_bioconductor-summarizedexperiment_r-base_r-data.table_pruned:531fb2eba7f91423' }"

    input:
    tuple val(meta), path(mae_res), path(gene_name_mapping)
    val(allelic_ratio_cutoff)
    val(padj_cutoff)
    val(max_cohort_freq)

    output:
    // results
    tuple val(meta), path("*.tsv.gz")                           , emit: res_all
    tuple val(meta), path("*${task.ext.prefix ?: meta.id}.tsv") , emit: res_signif
    tuple val(meta), path("*_rare.tsv")                         , emit: res_signif_rare
    // QC plots and summary
    tuple val(meta), path("mae_overview_mqc.tsv")               , emit: mae_overview
    tuple val(meta), path("cascade_plot_mqc.png")               , emit: cascade_plot
    tuple val(meta), path("variant_frequency_mqc.png")          , emit: variant_frequency
    tuple val(meta), path("median_of_each_category_mqc.tsv")    , emit: median_of_each_category
    tuple val(meta), path("results_mqc.tsv")                    , emit: results

    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Results.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    a <- file("MAE_results_all_${prefix}.tsv", "w")
    close(a)
    R.utils::gzip("MAE_results_all_${prefix}.tsv")
    file.remove("MAE_results_all_${prefix}.tsv")
    a <- file("MAE_results_${prefix}.tsv", "w")
    close(a)
    a <- file("MAE_results_${prefix}_rare.tsv", "w")
    close(a)
    a <- file("mae_overview_mqc.tsv", "w")
    close(a)
    a <- file("cascade_plot_mqc.png", "w")
    close(a)
    a <- file("variant_frequency_mqc.png", "w")
    close(a)
    a <- file("median_of_each_category_mqc.tsv", "w")
    close(a)
    a <- file("results_mqc.tsv", "w")
    close(a)
    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-r.utils:', as.character(packageVersion('R.utils'))),
            paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
            paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
            paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
            paste('    r-dt:', as.character(packageVersion('DT'))),
            paste('    bioconductor-genomicranges:', as.character(packageVersion('GenomicRanges'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment')))
        ),
    'versions.yml')
    """
}
