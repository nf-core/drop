process MAEQC_DNARNADESEQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a1fe185a170643337ce71d84c9262cf0fe44adb3f2fa51c5cd6007aebeebe18/data' :
        'community.wave.seqera.io/library/bioconductor-genomicranges_r-tmae_r-base_r-r.utils_r-stringr:e9d522a02847448b' }"

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
