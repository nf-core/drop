include { COUNTREADS } from '../../../modules/local/countreads/'

workflow ABERRANTEXPRESSION {
    take:
    bams            // queue channel: [ val(meta), path(bam), path(bai) ]
    count_ranges    // queue channel: [ val(meta), path(count_ranges) ]

    main:
    def ch_versions = Channel.empty()
    // TODO write the subworkflow

    def countreads_input = bams.map { meta, bam, bai ->
            [ [id: meta.gene_annotation], meta, bam, bai ]
        }
        .combine(count_ranges, by: 0)
        .map { _gene_meta, meta, bam, bai, count_ranges_ ->
            [ meta, bam, bai, meta.strand, meta.count_mode, meta.paired_end, meta.count_overlaps, count_ranges_ ]
        }

    COUNTREADS(
        countreads_input,
        params.ae_yield_size
    )
    ch_versions = ch_versions.mix(COUNTREADS.out.versions.first())

    COUNTREADS.out.counts.view()

    emit:
    versions = ch_versions
}
