process MAE_ALLELICCOUNTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/df/df2fe51e455c97d1598ba04a7d0cf8805a51eee5debda3e608f3a64ff72a261f/data' :
        'community.wave.seqera.io/library/bcftools_gatk4_htslib_samtools:949037bf07506ba9' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai), val(id), path(fasta), path(fai), path(dict)
    val(ignore_header_check)
    path(ncbi2ucsc)
    path(ucsc2ncbi)

    output:
    tuple val(meta), path("*.csv.gz")   , emit: csv
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Set temp dir in work directory
    mkdir ./tmp
    export TMPDIR=./tmp

    ASEReadCounter.sh \\
        $ncbi2ucsc \\
        $ucsc2ncbi \\
        $vcf \\
        $bam \\
        $id \\
        $fasta \\
        $ignore_header_check \\
        ${prefix}.csv.gz \\
        bcftools \\
        samtools \\
        gatk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
