include { COUNTREADS } from '../../../modules/local/countreads/'

workflow ABERRANTEXPRESSION {
    take:
    bams            // queue channel: [ val(meta), path(bam), path(bai) ]
    count_ranges    // queue channel: [ val(meta), path(count_ranges) ]

    include_groups  // list:          A list of groups to exclude from the aberrant expression analysis

    main:
    def ch_versions = Channel.empty()

    def bams_to_analyse = Channel.empty()
    def exclude_ids = Channel.value([])
    if (include_groups) {
        // Filter out the BAM files that don't have a group in the include_groups list
        def include_branch = bams.branch { meta, _bam, _bai ->
            yes: meta.drop_group.tokenize(",").intersect(include_groups).size() > 0
            no: true
        }
        // Get the IDs of the BAM files that are excluded
        exclude_ids = include_branch.no.collect { meta, _bam, _bai ->
            meta.id
        }
        // Get the BAM files that are included
        bams_to_analyse = include_branch.yes
    } else {
        bams_to_analyse = bams
    }

    def countreads_input = bams_to_analyse.map { meta, bam, bai ->
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

