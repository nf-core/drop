process OUTRIDER_PVALS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7ce5b3968fffff40e5f0ee411b06d26c117985c558d976181e5b3d8575bf3b89/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_bioconductor-outrider_bioconductor-summarizedexperiment_pruned:536d062a9b8ef2c5' }"

    input:
    tuple val(meta) , path(ods_fitted), val(ids), path(gene_name_mapping)
    tuple val(meta3), path(genes_to_test)
    val(implementation)
    path(parse_subsets_for_FDR) // Pass ${projectDir}/assets/helpers/parse_subsets_for_FDR.R to this input

    output:
    tuple val(meta), path("*.Rds")  , emit: ods_with_pvals
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'pvalsOutrider.R'

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
            paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    r-ggplot2:', as.character(packageVersion('ggplot2'))),
            paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
            paste('    r-tools:', as.character(packageVersion('tools'))),
            paste('    r-yaml:', as.character(packageVersion('yaml'))),
            paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment')))
        ),
    'versions.yml')
    """
}
