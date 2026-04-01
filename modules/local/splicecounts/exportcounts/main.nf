process SPLICECOUNTS_EXPORTCOUNTS {
    tag "${meta.id}"
    label 'process_single'

    // TODO discuss if this container should be split up, library loading should be done in the module script instead in that case
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/4319818c9416d02ea4e4bd81bad32533ce5268b81a5c53b7117984edad6655f1/data' :
        'community.wave.seqera.io/library/bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-bsgenome_bioconductor-delayedmatrixstats_pruned:e8af96bc80316fdc' }"

    input:
    tuple val(meta), path(annotation), path(splice_metrics, stageAs: "savedObjects/*"), val(drop_group)
    val(fraser_version)

    output:
    tuple val(meta), path("k_*_counts.tsv.gz") , emit: k_counts
    tuple val(meta), path("n_*_counts.tsv.gz") , emit: n_counts
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'exportCounts.R'

    stub:
    """
    #!/usr/bin/env Rscript
    a <- gzfile("k_j_counts.tsv.gz", "w")
    close(a)
    a <- gzfile("k_theta_counts.tsv.gz", "w")
    close(a)
    a <- gzfile("n_psi5_counts.tsv.gz", "w")
    close(a)
    a <- gzfile("n_psi3_counts.tsv.gz", "w")
    close(a)
    a <- gzfile("n_theta_counts.tsv.gz", "w")
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
            paste('    bioconductor-annotationdbi:', as.character(packageVersion('AnnotationDbi'))),
            paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
            paste('    bioconductor-delayedmatrixstats:', as.character(packageVersion('DelayedMatrixStats'))),
            paste('    bioconductor-bsgenome:', as.character(packageVersion('BSgenome'))),
            paste('    bioconductor-fraser:', as.character(packageVersion('FRASER')))
        ),
    'versions.yml')
    """
}
