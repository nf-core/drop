process MAE_DESEQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ params.add_af ? 'nf-core/deseqmae:1.0.4' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/4561f708be5a850be48042f0d4d21b2a121fa77576604832da55f650da5735b5/data' :
        'community.wave.seqera.io/library/r-tmae_r-base_r-r.utils_r-stringr:f69c262e49a5a245' }"

    input:
    tuple val(meta), path(counts)
    val(add_af)
    val(assembly)
    val(max_af)

    output:
    tuple val(meta), path("*.Rds")  , emit: res
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq_mae.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    a <- file("${prefix}_res.Rds", "w")
    close(a)

    ## VERSIONS FILE
    # Reinclude all versions in the stub when the tests can be run with the original container
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3])
            # paste('    r-r.utils:', as.character(packageVersion('R.utils'))),
            # paste('    r-string:', as.character(packageVersion('stringr'))),
            # paste('    r-tmae:', as.character(packageVersion('tMAE')))
        ),
    'versions.yml')
    """
}
