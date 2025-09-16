process GENECOUNTS_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7ce5b3968fffff40e5f0ee411b06d26c117985c558d976181e5b3d8575bf3b89/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_bioconductor-outrider_bioconductor-summarizedexperiment_pruned:536d062a9b8ef2c5' }"

    input:
    tuple val(meta), path(ods), path(bam_cov)

    output:
    tuple val(meta), path("sample_count_mqc.tsv")       , emit: sample_count
    tuple val(meta), path("read_counts_mqc.png")        , emit: read_counts           , optional: true
    tuple val(meta), path("mapped_vs_counted_mqc.png")  , emit: mapped_vs_counted
    tuple val(meta), path("size_factors_mqc.png")       , emit: size_factors
    tuple val(meta), path("reads_statistics_mqc.tsv")   , emit: reads_statistics
    tuple val(meta), path("meanCounts_mqc.png")         , emit: meanCounts
    tuple val(meta), path("expressedGenes_mqc.png")     , emit: expressedGenes
    tuple val(meta), path("expressed_genes_mqc.tsv")    , emit: expressed_genes
    tuple val(meta), path("expression_sex_mqc.tsv")     , emit: expression_sex
    tuple val(meta), path("sex_matched_mqc.png")        , emit: sex_matched           , optional: true
    tuple val(meta), path("xist_uty.tsv")               , emit: xist_uty              , optional: true
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Summary.R'

    stub:
    """
    #!/usr/bin/env Rscript
    a <- file("sample_count_mqc.tsv", "w")
    close(a)
    a <- file("mapped_vs_counted_mqc.png", "w")
    close(a)
    a <- file("size_factors_mqc.png", "w")
    close(a)
    a <- file("reads_statistics_mqc.tsv", "w")
    close(a)
    a <- file("meanCounts_mqc.png", "w")
    close(a)
    a <- file("expressedGenes_mqc.png", "w")
    close(a)
    a <- file("expressed_genes_mqc.tsv", "w")
    close(a)
    a <- file("expression_sex_mqc.tsv", "w")
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
            paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
            paste('    bioconductor-GenomicAlignments::', as.character(packageVersion('GenomicAlignments'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment'))),
            paste('    bioconductor-outrider:', as.character(packageVersion('OUTRIDER')))
        ),
    'versions.yml')
    """
}
