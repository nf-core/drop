process MAE_DESEQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'nf-core/deseqmae:1.0.4'

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

    a <- file("${prefix}.Rds", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-r.utils:', as.character(packageVersion('R.utils'))),
            paste('    r-string:', as.character(packageVersion('stringr'))),
            paste('    r-tmae:', as.character(packageVersion('tMAE'))),
            paste('    bioconductor-mafdb.gnomad-r2.1.grch38:', as.character(packageVersion('MafDb.gnomAD.r2.1.GRCh38'))),
            paste('    bioconductor-mafdb.gnomad-r2.1.hs37d5:', as.character(packageVersion('MafDb.gnomAD.r2.1.hs37d5')))
        ),
    'versions.yml')
    """
}
