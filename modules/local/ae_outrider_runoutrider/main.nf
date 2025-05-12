process AE_OUTRIDER_RUNOUTRIDER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3c/3c54f46cc4a1f9f92c216095525003cce5ad6581696086a76d9df45767d45a48/data' :
        'community.wave.seqera.io/library/bioconductor-genomicfeatures_bioconductor-rtracklayer_bioconductor-txdbmaker_r-base_pruned:2d189ec410925335' }"

    input:
    tuple val(meta), path(ods)

    output:
    tuple val(meta), path("*ods_fitted.Rds")            , emit: ods_fitted
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'runOutrider.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(GenomicFeatures)
        library(data.table)
        library(rtracklayer)
        library(magrittr)
    })
  
    """
}
