process SPLICECOUNTS_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/4319818c9416d02ea4e4bd81bad32533ce5268b81a5c53b7117984edad6655f1/data' :
        'community.wave.seqera.io/library/bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-bsgenome_bioconductor-delayedmatrixstats_pruned:e8af96bc80316fdc' }"

    input:
    tuple val(meta), path(fdsobj_raw_local, stageAs: "savedObjects/*"), path(fdsobj_raw, stageAs: "savedObjects/*"), val(drop_group)
    val(fraser_version)
    path(config)

    output:
    tuple val(meta), path("number_of_samples_mqc.tsv")              , emit: number_of_samples
    tuple val(meta), path("number_of_introns_mqc.tsv")              , emit: number_of_introns
    tuple val(meta), path("number_of_splice_sites_mqc.tsv")         , emit: number_of_splice_sites
    tuple val(meta), path("comparison_local_and_external_mqc.png")  , emit: comparison_local_and_external , optional: true
    tuple val(meta), path("expression_filtering_mqc.png")           , emit: expression_filtering
    tuple val(meta), path("variability_filtering_mqc.png")          , emit: variability_filtering
    path  "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'Summary.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- file("number_of_samples_mqc.tsv", "w")
    close(a)
    a <- file("number_of_introns_mqc.tsv", "w")
    close(a)
    a <- file("number_of_splice_sites_mqc.tsv", "w")
    close(a)
    a <- file("expression_filtering_mqc.png", "w")
    close(a)
    a <- file("variability_filtering_mqc.png", "w")
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
            paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
            paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
            paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
            paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
        ),
    'versions.yml')
    """
}
