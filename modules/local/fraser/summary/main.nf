process FRASER_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/4319818c9416d02ea4e4bd81bad32533ce5268b81a5c53b7117984edad6655f1/data' :
        'community.wave.seqera.io/library/bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-bsgenome_bioconductor-delayedmatrixstats_pruned:e8af96bc80316fdc' }"

    input:
    tuple val(meta), path(fdsin, stageAs: 'input/savedObjects/*'), path(results), val(drop_group), val(annotation_id)
    val(padj_cutoff)
    val(deltaPsi_cutoff)
    val(fraser_version)
    path(config)

    output:
    tuple val(meta), path("fraser_overview_mqc.tsv")           , emit: fraser_overview
    tuple val(meta), path("q_estimation*mqc.png")              , emit: q_estimation
    tuple val(meta), path("aberrantly_spliced_genes_mqc.png")  , emit: aberrantly_spliced_genes
    tuple val(meta), path("batch_correlation*mqc.png")         , emit: batch_correlation
    tuple val(meta), path("total_outliers_mqc.tsv")            , emit: total_outliers
    tuple val(meta), path("results_mqc.tsv")                   , emit: results
    path  "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Summary.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- file("fraser_overview_mqc.tsv", "w")
    close(a)
    a <- file("q_estimation_x_mqc.png", "w")
    close(a)
    a <- file("aberrantly_spliced_genes_mqc.png", "w")
    close(a)
    a <- file("batch_correlation_x_mqc.png", "w")
    close(a)
    a <- file("total_outliers_mqc.tsv", "w")
    close(a)
    a <- file("results_mqc.tsv", "w")
    close(a)

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
            paste('    r-rcolorbrewer:', as.character(packageVersion('RColorBrewer'))),
            paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
            paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
            paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
            paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
        ),
    'versions.yml')
    """
}
