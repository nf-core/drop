include { COUNTREADS                         } from '../../../modules/local/countreads/'
include { MERGECOUNTS                        } from '../../../modules/local/mergecounts/'
include { FILTERCOUNTS                       } from '../../../modules/local/filtercounts/'
include { COUNTEXPRESSION_SUMMARY            } from '../../../modules/local/countexpression/summary/'
include { OUTRIDER_RUN                       } from '../../../modules/local/outrider/run'
include { OUTRIDER_PVALS                     } from '../../../modules/local/outrider/pvals'
include { OUTRIDER_RESULTS                   } from '../../../modules/local/outrider/results'
include { OUTRIDER_SUMMARY                   } from '../../../modules/local/outrider/summary'
include { MULTIQC as MULTIQC_COUNTEXPRESSION } from '../../../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_OUTRIDER        } from '../../../modules/nf-core/multiqc/main'

include { BAM_STATS_IDXSTATS_MERGE           } from "../bam_stats_idxstats_merge/"

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

    def external_counts_ids = inputs_branch.counts_present
        .map { meta, _gene_counts ->
            [[meta.id, meta.drop_group]]
        }
        .collect()
        .map { list ->
            def group_ids = [:]
            list.each { sample ->
                sample[1].tokenize(",").each { group ->
                    group_ids[group] = group_ids.get(group, []) + sample[0]
                }
            }
            return group_ids
        }
        .ifEmpty([:])

    def countreads_input = inputs_branch.counts_missing
        .combine(count_ranges)
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
        .combine(external_counts_ids)
        .map { _gene_meta, meta, counts, count_ranges_, external_counts_ids_ ->
            def unique_counts = counts.unique() // Remove duplicates in case a merged count file was given as an external count for multiple samples
            [ meta, unique_counts, count_ranges_, external_counts_ids_.get(meta.drop_group, []) ]
        }

    MERGECOUNTS(
        mergereads_input,
        samplesheet
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

    //
    // summary for expression counts
    //
    def countexpression_summary_input = FILTERCOUNTS.out.output
                    .map { meta, ods_unfitted ->
                        [ meta.drop_group,      // key for join
                        meta,
                        ods_unfitted ]
                    }
                    .combine(BAM_STATS_IDXSTATS_MERGE.out.merged_bam_stats, by: 0)
                    .map { group, meta, ods_unfitted, bamstats ->
                        [ meta, ods_unfitted, bamstats ]
                    }

    COUNTEXPRESSION_SUMMARY(
        countexpression_summary_input
    )
    ch_versions = ch_versions.mix(COUNTEXPRESSION_SUMMARY.out.versions.first())

    //
    // summary for outrider results
    //
    OUTRIDER_PVALS.out.ods_with_pvals.view { println "[DEBUG OUTRIDER_PVALS] $it" }
    OUTRIDER_RESULTS.out.results.view { println "[DEBUG OUTRIDER_RESULTS] $it" }

    def outrider_summary_input = OUTRIDER_PVALS.out.ods_with_pvals
                    .join(OUTRIDER_RESULTS.out.results, by: 0)
                    .map { meta, ods, results ->
                            [ meta, ods, results, meta.drop_group, meta.id]
                        }
    outrider_summary_input.view { println "[DEBUG OUTRIDER_SUMMARY] $it" }

    OUTRIDER_SUMMARY(
        outrider_summary_input,
        params.ae_padj_cutoff,
        params.ae_z_score_cutoff
    )
    ch_versions = ch_versions.mix(COUNTEXPRESSION_SUMMARY.out.versions.first())

    //
    // multiqc for counting
    //
    def multiqc_countexpression_input = COUNTEXPRESSION_SUMMARY.out.sample_count
        .join(COUNTEXPRESSION_SUMMARY.out.read_counts, remainder: true)
        .join(COUNTEXPRESSION_SUMMARY.out.mapped_vs_counted)
        .join(COUNTEXPRESSION_SUMMARY.out.size_factors)
        .join(COUNTEXPRESSION_SUMMARY.out.reads_statistics)
        .join(COUNTEXPRESSION_SUMMARY.out.meanCounts)
        .join(COUNTEXPRESSION_SUMMARY.out.expressedGenes)
        .join(COUNTEXPRESSION_SUMMARY.out.expressed_genes)
        .join(COUNTEXPRESSION_SUMMARY.out.expression_sex, remainder: true)
        .join(COUNTEXPRESSION_SUMMARY.out.sex_matched, remainder: true)
        .map { it.drop(1).findAll { p -> p instanceof Path } }   // drop meta and null

    MULTIQC_COUNTEXPRESSION(
        multiqc_countexpression_input,
        Channel.fromPath("$projectDir/assets/multiqc_countexpression_config.yml", checkIfExists: true).collect(),
        [],
        [],
        [],
        []
    )

    //
    // multiqc for outrider
    //
    def multiqc_outrider_input = OUTRIDER_SUMMARY.out.outrider_overview
        .join(OUTRIDER_SUMMARY.out.encoding_dimension)
        .join(OUTRIDER_SUMMARY.out.aberrant_genes_per_sample)
        .join(OUTRIDER_SUMMARY.out.batch_correction_raw)
        .join(OUTRIDER_SUMMARY.out.batch_correction_normalized)
        .join(OUTRIDER_SUMMARY.out.geneSampleHeatmap_raw)
        .join(OUTRIDER_SUMMARY.out.geneSampleHeatmap_normalized)
        .join(OUTRIDER_SUMMARY.out.bcv)
        .join(OUTRIDER_SUMMARY.out.outrider_result_overview)
        .join(OUTRIDER_SUMMARY.out.aberrant_samples)
        .join(OUTRIDER_SUMMARY.out.significant_results)
        .map { it.drop(1).findAll { p -> p instanceof Path } }   // drop meta

    MULTIQC_COUNTEXPRESSION(
        multiqc_outrider_input,
        Channel.fromPath("$projectDir/assets/multiqc_outrider_config.yml", checkIfExists: true).collect(),
        [],
        [],
        [],
        []
    )


    emit:
    versions  = ch_versions
    results   = OUTRIDER_RESULTS.out.results
    count_report = MULTIQC_COUNTEXPRESSION.out.report.toList()
}
