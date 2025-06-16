include { COUNTREADS                } from '../../../modules/local/countreads/'
include { MERGECOUNTS               } from '../../../modules/local/mergecounts/'
include { FILTERCOUNTS              } from '../../../modules/local/filtercounts/'
include { OUTRIDER_RUN              } from '../../../modules/local/outrider/run'
include { OUTRIDER_PVALS            } from '../../../modules/local/outrider/pvals'
include { OUTRIDER_RESULTS          } from '../../../modules/local/outrider/results'

include { BAM_STATS_IDXSTATS_MERGE  } from "../bam_stats_idxstats_merge/"

workflow ABERRANTEXPRESSION {
    take:
    input               // queue channel: [ val(meta), path(bam), path(bai), path(gene_counts) ]
    count_ranges        // queue channel: [ val(meta), path(count_ranges) ]
    txdb                // queue channel: [ val(meta), path(txdb) ]
    gene_name_mapping   // queue channel: [ val(meta), path(gene_name_mapping) ]
    samplesheet         // value channel: [ val(meta), path(samplesheet) ]
    genes_to_test       // value channel: [ val(meta), path(genes_to_test) ]
    hpo                 // value channel: [ val(meta), path(hpo) ]

    include_groups      // list:          A list of groups to include in the aberrant expression analysis

    main:
    def ch_versions = Channel.empty()

    def inputs_to_analyse = Channel.empty()
    if (include_groups) {
        // Filter out the BAM files that don't have a group in the include_groups list
        def include_branch = input.branch { meta, _bam, _bai, _gene_counts ->
            yes: meta.drop_group.tokenize(",").intersect(include_groups).size() > 0
            no: true
        }

        // Get the BAM files that are included
        inputs_to_analyse = include_branch.yes
    } else {
        inputs_to_analyse = input
    }

    //
    // Count the reads in the BAM files
    //

    def inputs_branch = inputs_to_analyse.branch { meta, bam, bai, gene_counts ->
        counts_present: gene_counts
            return [ meta, gene_counts ]
        counts_missing: true
            return [ meta, bam, bai ]
    }

    def external_counts_ids = inputs_branch.counts_present.map { meta, _gene_counts ->
        meta.id
    }.collect().ifEmpty([])

    def countreads_input = inputs_branch.counts_missing
        .combine(count_ranges)
        .filter { meta, _bam, _bai, annotation_meta, _count_ranges ->
            // Determine which gene annotation versions should be used for each sample
            meta.gene_annotation == "" || meta.gene_annotation == annotation_meta.id
        }
        .map { meta, bam, bai, annotation_meta, count_ranges_ ->
            def new_meta = meta + [gene_annotation: annotation_meta.id]
            [ new_meta, bam, bai, meta.strand, meta.count_mode, meta.paired_end, meta.count_overlaps, count_ranges_ ]
        }

    COUNTREADS(
        countreads_input,
        params.ae_yield_size
    )
    ch_versions = ch_versions.mix(COUNTREADS.out.versions.first())

    //
    // Merge the counts per gene annotation and drop group
    //

    // Group samples based on gene annotation and drop group
    def mergereads_input = COUNTREADS.out.counts
        .mix(inputs_branch.counts_present)
        .map { meta, count ->
            [ meta, count, meta.drop_group.tokenize(",").intersect(include_groups) ]
        }
        .transpose(by: 2) // Split the entries per included group
        .map { meta, count, group ->
            def new_meta = [ id:meta.gene_annotation, drop_group:group ]
            [ groupKey(new_meta, meta.drop_group_ann_counts.get(group).get(meta.gene_annotation)), count, meta.id ]
        }
        .groupTuple() // Group the counts by gene annotation and drop group
        .map { meta, counts, ids ->
            def new_meta = meta + [ ids:ids ]
            [ [id: meta.id], new_meta, counts ]
        }
        .combine(count_ranges, by: 0)
        .map { _gene_meta, meta, counts, count_ranges_ ->
            [ meta, counts, count_ranges_ ]
        }

    MERGECOUNTS(
        mergereads_input,
        samplesheet,
        external_counts_ids
    )
    ch_versions = ch_versions.mix(MERGECOUNTS.out.versions.first())

    //
    // Filter the merged counts
    //

    def filtercounts_input = MERGECOUNTS.out.output
        .map { meta, counts ->
            [ [ id: meta.id ], meta, counts ]
        }
        .combine(txdb, by: 0)
        .map { _gene_meta, meta, counts, txdb_ ->
            [ meta, counts, txdb_ ]
        }

    FILTERCOUNTS(
        filtercounts_input,
        params.ae_fpkm_cutoff
    )
    ch_versions = ch_versions.mix(FILTERCOUNTS.out.versions.first())

    //
    // Run OUTRIDER
    //

    OUTRIDER_RUN(
        FILTERCOUNTS.out.output,
        params.ae_implementation,
        params.ae_max_tested_dimension_proportion
    )
    ch_versions = ch_versions.mix(OUTRIDER_RUN.out.versions.first())

    //
    // Calculate p-values of the OUTRIDER output
    //

    def outrider_pvals_input = OUTRIDER_RUN.out.ods_fitted
        .map { meta, ods_fitted ->
            [ [id:meta.id], meta, ods_fitted ]
        }
        .combine(gene_name_mapping, by:0)
        .map { _gene_meta, meta, ods_fitted, gene_name_mapping_ ->
            [ meta, ods_fitted, meta.ids, gene_name_mapping_ ]
        }

    OUTRIDER_PVALS(
        outrider_pvals_input,
        genes_to_test,
        params.ae_implementation,
        file("${projectDir}/assets/helpers/parse_subsets_for_FDR.R")
    )
    ch_versions = ch_versions.mix(OUTRIDER_PVALS.out.versions.first())

    //
    // Get the OUTRIDER results
    //

    def outrider_results_input = OUTRIDER_PVALS.out.ods_with_pvals
        .map { meta, ods_with_pvals ->
            [ [id:meta.id], meta, ods_with_pvals ]
        }
        .combine(gene_name_mapping, by:0)
        .map { _gene_meta, meta, ods_with_pvals, gene_name_mapping_ ->
            [ meta, ods_with_pvals, gene_name_mapping_ ]
        }

    OUTRIDER_RESULTS(
        outrider_results_input,
        samplesheet,
        hpo,
        params.ae_padj_cutoff,
        params.ae_z_score_cutoff,
        file("${projectDir}/assets/helpers/add_HPO_cols.R")
    )
    ch_versions = ch_versions.mix(OUTRIDER_RESULTS.out.versions.first())

    //
    // Calculate and merge the BAM stats
    //

    BAM_STATS_IDXSTATS_MERGE(
        inputs_to_analyse.map { meta, bam, bai, _gene_counts ->
            [ meta, bam, bai ]
        }.filter { _meta, bam, _bai -> bam },
        include_groups
    )
    ch_versions = ch_versions.mix(BAM_STATS_IDXSTATS_MERGE.out.versions)

    emit:
    versions  = ch_versions
    results   = OUTRIDER_RESULTS.out.results
    bam_stats = BAM_STATS_IDXSTATS_MERGE.out.merged_bam_stats
}

