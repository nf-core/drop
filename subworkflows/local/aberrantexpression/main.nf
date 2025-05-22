include { COUNTREADS   } from '../../../modules/local/countreads/'
include { MERGECOUNTS  } from '../../../modules/local/mergecounts/'
include { FILTERCOUNTS } from '../../../modules/local/filtercounts/'
include { OUTRIDER_RUN } from '../../../modules/local/outrider/run'

workflow ABERRANTEXPRESSION {
    take:
    bams            // queue channel: [ val(meta), path(bam), path(bai) ]
    count_ranges    // queue channel: [ val(meta), path(count_ranges) ]
    txdb            // queue channel: [ val(meta), path(txdb) ]
    samplesheet     // value channel: [ val(meta), path(samplesheet) ]

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
        }.ifEmpty([])
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

    // Group samples based on gene annotation
    def mergereads_input = COUNTREADS.out.counts
        .map { meta, count ->
            [ [id: meta.gene_annotation], count ]
        }
        .groupTuple()
        .combine(count_ranges, by: 0)

    MERGECOUNTS(
        mergereads_input,
        samplesheet,
        exclude_ids
    )
    ch_versions = ch_versions.mix(MERGECOUNTS.out.versions.first())

    def filtercounts_input = MERGECOUNTS.out.output
        .combine(txdb, by: 0)

    FILTERCOUNTS(
        filtercounts_input,
        params.ae_fpkm_cutoff
    )
    ch_versions = ch_versions.mix(FILTERCOUNTS.out.versions.first())

    FILTERCOUNTS.out.output.view()

    OUTRIDER_RUN(
        FILTERCOUNTS.out.output,
        params.ae_implementation,
        params.ae_max_tested_dimension_proportion
    )
    ch_versions = ch_versions.mix(OUTRIDER_RUN.out.versions.first())

    OUTRIDER_RUN.out.ods_fitted.view()

    emit:
    versions = ch_versions
}

