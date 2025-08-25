process SPLICECOUNTS_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c87aab5452ff7d6a9344148e284ca01a8c5d2239d7b0f26164da1806a8c59875/data' :
        'community.wave.seqera.io/library/bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-bsgenome_bioconductor-delayedmatrixstats_pruned:b2aa4004c588aaaf' }"

    input:
    tuple val(meta), path(fdsobj_raw_local, stageAs: "savedObjects/*"), path(fdsobj_raw, stageAs: "savedObjects/*"), val(drop_group)
    val(fraser_version)
    path(config)

    output:
    // tuple val(meta), path("sample_count_mqc.tsv")       , emit: sample_count
    // tuple val(meta), path("read_counts_mqc.png")        , emit: read_counts           , optional: true
    // tuple val(meta), path("mapped_vs_counted_mqc.png")  , emit: mapped_vs_counted
    // tuple val(meta), path("size_factors_mqc.png")       , emit: size_factors
    // tuple val(meta), path("reads_statistics_mqc.tsv")   , emit: reads_statistics
    // tuple val(meta), path("meanCounts_mqc.png")         , emit: meanCounts
    // tuple val(meta), path("expressedGenes_mqc.png")     , emit: expressedGenes
    // tuple val(meta), path("expressed_genes_mqc.tsv")    , emit: expressed_genes
    // tuple val(meta), path("expression_sex_mqc.tsv")     , emit: expression_sex
    // tuple val(meta), path("sex_matched_mqc.png")        , emit: sex_matched           , optional: true
    // tuple val(meta), path("xist_uty.tsv")               , emit: xist_uty              , optional: true
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Summary.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-rmarkdown:', as.character(packageVersion('rmarkdown'))),
            paste('    r-knitr:', as.character(packageVersion('knitr'))),
            paste('    r-devtools:', as.character(packageVersion('devtools'))),
            paste('    r-yaml:', as.character(packageVersion('yaml'))),
            paste('    r-bbmisc:', as.character(packageVersion('BBmisc'))),
            paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
            paste('    r-plotly:', as.character(packageVersion('plotly'))),
            paste('    r-rhdf5:', as.character(packageVersion('rhdf5'))),
            paste('    r-cowplot:', as.character(packageVersion('cowplot'))),
            paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
            paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
            paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
            paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
        ),
    'versions.yml')
    """
}
