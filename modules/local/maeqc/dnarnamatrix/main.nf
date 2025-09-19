process MAEQC_DNARNAMATRIX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/40991c6312e0419255336842c526cdcc20fa6cf3a5475749ad579a97f5173220/data' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-genomicranges_bioconductor-variantannotation_r-base_pruned:3528952457d92839' }"

    input:
    tuple val(meta) , path(res), path(vcfs), path(tbis), path(samplesheet), val(drop_group)
    tuple val(meta2), path(qc_vcf), path(qc_tbi)

    output:
    tuple val(meta), path("*.Rds")  , emit: mat_qc
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'create_matrix_dna_rna_cor.R'

    stub:
    """
    #!/usr/bin/env Rscript
    a <- file("dna_rna_qc_matrix.Rds", "w")
    close(a)

    ## VERSIONS FILE
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
            paste('    r-magrittr:', as.character(packageVersion('magrittr'))),
            paste('    r-data.table:', as.character(packageVersion('data.table'))),
            paste('    bioconductor-biocparallel:', as.character(packageVersion('BiocParallel'))),
            paste('    bioconductor-variantannotation:', as.character(packageVersion('VariantAnnotation'))),
            paste('    bioconductor-genomicranges:', as.character(packageVersion('GenomicRanges')))
        ),
    'versions.yml')
    """
}
