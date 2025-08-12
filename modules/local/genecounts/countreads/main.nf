process GENECOUNTS_COUNTREADS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba7130887847376e2da9c9325c3c9c5531e62c9811b2467567ea78bfeaf9c73c/data' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-genomicalignments_bioconductor-rsamtools_r-base_r-data.table:195b16b35d0c9f89' }"

    input:
    tuple val(meta), path(bam), path(bai), val(strand), val(count_mode), val(paired_end), val(count_overlaps), path(count_ranges)
    val(yield_size)

    output:
    tuple val(meta), path("*.Rds")  , emit: counts
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'countReads.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- file("${prefix}.Rds", "w")
    close(a)

    # Write the versions file
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    bioconductor-genomicalignments:', as.character(packageVersion('GenomicAlignments'))),
            paste('    bioconductor-rsamtools:', as.character(packageVersion('Rsamtools'))),
            paste('    bioconductor-biocparallel:', as.character(packageVersion('BiocParallel')))
        ),
    'versions.yml')
    """
}
