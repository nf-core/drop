include { GENECOUNTS_COUNTREADS              } from '../../../modules/local/genecounts/countreads/'
include { GENECOUNTS_MERGECOUNTS             } from '../../../modules/local/genecounts/mergecounts/'
include { GENECOUNTS_FILTERCOUNTS            } from '../../../modules/local/genecounts/filtercounts/'
include { GENECOUNTS_SUMMARY                 } from '../../../modules/local/genecounts/summary/'
include { OUTRIDER_FIT                       } from '../../../modules/local/outrider/fit'
include { OUTRIDER_PVALS                     } from '../../../modules/local/outrider/pvals'
include { OUTRIDER_RESULTS                   } from '../../../modules/local/outrider/results'
include { OUTRIDER_SUMMARY                   } from '../../../modules/local/outrider/summary'
include { MULTIQC as MULTIQC_GENECOUNTS      } from '../../../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_OUTRIDER        } from '../../../modules/nf-core/multiqc/main'
include { BAM_STATS_IDXSTATS_MERGE           } from "../bam_stats_idxstats_merge/"
include { STRIP_NAV_MULTIQC as STRIP_NAV_GENECOUNTS     } from '../../../modules/local/strip_nav_multiqc/'
include { STRIP_NAV_MULTIQC as STRIP_NAV_OUTRIDER       } from '../../../modules/local/strip_nav_multiqc/'

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

    def GENECOUNTS_COUNTREADS_input = inputs_branch.counts_missing
        .combine(count_ranges)
        .map { meta, bam, bai, annotation_meta, count_ranges_ ->
            def new_meta = meta + [gene_annotation: annotation_meta.id]
            [ new_meta, bam, bai, meta.strand, meta.count_mode, meta.paired_end, meta.count_overlaps, count_ranges_ ]
        }

    GENECOUNTS_COUNTREADS(
        GENECOUNTS_COUNTREADS_input,
        params.ae_yield_size
    )
    ch_versions = ch_versions.mix(GENECOUNTS_COUNTREADS.out.versions.first())

    //
    // Merge the counts per gene annotation and drop group
    //

    // Group samples based on gene annotation and drop group
    def mergereads_input = GENECOUNTS_COUNTREADS.out.counts
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

    GENECOUNTS_MERGECOUNTS(
        mergereads_input,
        samplesheet
    )
    ch_versions = ch_versions.mix(GENECOUNTS_MERGECOUNTS.out.versions.first())

    //
    // Filter the merged counts
    //

    def filtercounts_input = GENECOUNTS_MERGECOUNTS.out.output
        .map { meta, counts ->
            [ [ id: meta.id ], meta, counts ]
        }
        .combine(txdb, by: 0)
        .map { _gene_meta, meta, counts, txdb_ ->
            [ meta, counts, txdb_ ]
        }

    GENECOUNTS_FILTERCOUNTS(
        filtercounts_input,
        params.ae_fpkm_cutoff
    )
    ch_versions = ch_versions.mix(GENECOUNTS_FILTERCOUNTS.out.versions.first())

    //
    // Run OUTRIDER
    //

    OUTRIDER_FIT(
        GENECOUNTS_FILTERCOUNTS.out.output,
        params.ae_implementation,
        params.ae_max_tested_dimension_proportion,
        params.random_seed,
        params.ae_use_grid_search_to_obtain_q
    )
    ch_versions = ch_versions.mix(OUTRIDER_FIT.out.versions.first())

    //
    // Calculate p-values of the OUTRIDER output
    //

    def outrider_pvals_input = OUTRIDER_FIT.out.ods_fitted
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
    def external_samples = inputs_to_analyse
        .filter { _meta, bam, _bai, gene_counts ->
            !bam && gene_counts                        // no BAM but a counts file
        }
        .map    { meta, _bam, _bai, _gene_counts ->
            meta
        }

    BAM_STATS_IDXSTATS_MERGE(
        inputs_to_analyse.map { meta, bam, bai, _gene_counts ->
                [ meta, bam, bai ]
            }.filter { _meta, bam, _bai -> bam },       // sample with bam files
        external_samples,                               // sample external counts
        include_groups
    )
    ch_versions = ch_versions.mix(BAM_STATS_IDXSTATS_MERGE.out.versions)

    //
    // summary for expression counts
    //
    def genecounts_summary_input = GENECOUNTS_FILTERCOUNTS.out.output
                    .map { meta, ods_unfitted ->
                        [ meta.drop_group,      // key for join
                        meta,
                        ods_unfitted ]
                    }
                    .combine(BAM_STATS_IDXSTATS_MERGE.out.merged_bam_stats, by: 0)
                    .map { _group, meta, ods_unfitted, bamstats ->
                        [ meta, ods_unfitted, bamstats ]
                    }

    GENECOUNTS_SUMMARY(
        genecounts_summary_input
    )
    ch_versions = ch_versions.mix(GENECOUNTS_SUMMARY.out.versions.first())

    //
    // summary for outrider results
    //
    def outrider_summary_input = OUTRIDER_PVALS.out.ods_with_pvals
                    .join(OUTRIDER_RESULTS.out.results, by: 0)
                    .map { meta, ods, results ->
                            [ meta, ods, results, meta.drop_group, meta.id]
                        }

    OUTRIDER_SUMMARY(
        outrider_summary_input,
        params.ae_padj_cutoff,
        params.ae_z_score_cutoff
    )
    ch_versions = ch_versions.mix(GENECOUNTS_SUMMARY.out.versions.first())

    //
    // multiqc for counting
    //
    def multiqc_genecounts_input = GENECOUNTS_SUMMARY.out.sample_count
        .join(GENECOUNTS_SUMMARY.out.read_counts, remainder: true)
        .join(GENECOUNTS_SUMMARY.out.mapped_vs_counted)
        .join(GENECOUNTS_SUMMARY.out.size_factors)
        .join(GENECOUNTS_SUMMARY.out.reads_statistics)
        .join(GENECOUNTS_SUMMARY.out.meanCounts)
        .join(GENECOUNTS_SUMMARY.out.expressedGenes)
        .join(GENECOUNTS_SUMMARY.out.expressed_genes)
        .join(GENECOUNTS_SUMMARY.out.expression_sex, remainder: true)
        .join(GENECOUNTS_SUMMARY.out.sex_matched, remainder: true)
        .map { tup ->
        def meta  = tup[0]
        def files = tup.drop(1).findAll { it instanceof Path } // drop meta and null
        def tag   = "${meta.id}__${meta.drop_group}" // tag for publishDir
        [ tag, files ]
    }
    .tap { genecounts_by_tag }

    def ch_extra_cfg_genecounts = genecounts_by_tag
        .collectFile { tag, _ ->
            def yaml = """
            output_fn_name: "[TAG:genecounts_${tag}]_multiqc_report.html"
            data_dir_name:  "[TAG:genecounts_${tag}]_multiqc_data"
            plots_dir_name: "[TAG:genecounts_${tag}]_multiqc_plots"
            """.stripIndent()
            [ "[TAG:genecounts_${tag}]_config.yml", yaml ]
        }
        .map { Path cfg ->
        def m = (cfg.getName() =~ /\[TAG:genecounts_([^\]]+)\]_config\.yml$/)
        def tag = m ? m[0][1] : cfg.baseName
        [ tag, cfg ]
    }

    def ch_genecounts_bundle = multiqc_genecounts_input
        .join(ch_extra_cfg_genecounts)
        .map { _tag, files, cfg -> [files, cfg] }

    MULTIQC_GENECOUNTS(
        ch_genecounts_bundle.map { it[0] },
        Channel.fromPath("$projectDir/assets/multiqc_configs/multiqc_genecounts_config.yml", checkIfExists: true).collect(),
        ch_genecounts_bundle.map { it[1] },
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
        .map { tup ->
        def meta  = tup[0]
        def files = tup.drop(1).findAll { it instanceof Path } // drop meta and null
        def tag   = "${meta.id}__${meta.drop_group}" // tag for publishDir
        [ tag, files ]
    }
    .tap { outrider_by_tag }

    def ch_extra_cfg_outrider = outrider_by_tag
        .collectFile { tag, _ ->
            def yaml = """
            output_fn_name: "[TAG:outrider_${tag}]_multiqc_report.html"
            data_dir_name:  "[TAG:outrider_${tag}]_multiqc_data"
            plots_dir_name: "[TAG:outrider_${tag}]_multiqc_plots"
            """.stripIndent()
            [ "[TAG:outrider_${tag}]_config.yml", yaml ]
        }
        .map { Path cfg ->
        def m = (cfg.getName() =~ /\[TAG:outrider_([^\]]+)\]_config\.yml$/)
        def tag = m ? m[0][1] : cfg.baseName
        [ tag, cfg ]
    }

    def ch_outrider_bundle = multiqc_outrider_input
        .join(ch_extra_cfg_outrider)
        .map { _tag, files, cfg -> [files, cfg] }

    MULTIQC_OUTRIDER(
        ch_outrider_bundle.map { it[0] },
        Channel.fromPath("$projectDir/assets/multiqc_configs/multiqc_outrider_config.yml", checkIfExists: true).collect(),
        ch_outrider_bundle.map { it[1] },
        [],
        [],
        []
    )

    STRIP_NAV_GENECOUNTS(MULTIQC_GENECOUNTS.out.report)
    STRIP_NAV_OUTRIDER(MULTIQC_OUTRIDER.out.report)

    emit:
    versions  = ch_versions
    results   = OUTRIDER_RESULTS.out.results
    count_report = STRIP_NAV_GENECOUNTS.out.report.toList()
    outrider_report = STRIP_NAV_OUTRIDER.out.report.toList()
}
