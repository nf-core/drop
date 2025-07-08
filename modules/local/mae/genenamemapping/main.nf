process MAE_GENENAMEMAPPING {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/82/82f7d0a8aa874088afd08f1f9f58761fbb5b0076626d6de714312dc43e902314/data' :
        'community.wave.seqera.io/library/bioconductor-rtracklayer_r-base_r-data.table_r-magrittr_r-tidyr:57d70ac44abd5d20' }"

    input:
    tuple val(meta) , path(gtf)

    output:
    tuple val(meta), path("*.tsv")  , emit: tsv
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'gene_name_mapping.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    a <- file("gene_name_mapping_${prefix}.tsv", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-tidyr:', as.character(packageVersion('tidyr'))),
            paste('    bioconductor-rtracklayer:', as.character(packageVersion('rtracklayer')))
        ),
    'versions.yml')
    """
}
