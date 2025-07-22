process MAEQC_DESEQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a2/a20e81a37e3e4376534cf7a225da72640b35674961fb3317e4e72c0177d0ac32/data' :
        'community.wave.seqera.io/library/bioconductor-genomicranges_r-tmae_r-base_r-r.utils:153fb7526c2a316a' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.Rds")  , emit: ods_with_pvals
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq_qc.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- file("${prefix}.Rds", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-tmae:', as.character(packageVersion('tMAE'))),
            paste('    bioconductor-genomicranges:', as.character(packageVersion('GenomicRanges')))
        ),
    'versions.yml')
    """
}
