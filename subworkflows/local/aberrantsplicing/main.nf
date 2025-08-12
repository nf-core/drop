include { SPLICECOUNTS_INIT                     } from '../../../modules/local/splicecounts/init'
include { SPLICECOUNTS_COUNTSPLITREADS          } from '../../../modules/local/splicecounts/countsplitreads'
include { SPLICECOUNTS_MERGESPLITREADS          } from '../../../modules/local/splicecounts/mergesplitreads'
include { SPLICECOUNTS_COUNTNONSPLITREADS       } from '../../../modules/local/splicecounts/countnonsplitreads'
include { SPLICECOUNTS_MERGENONSPLITREADS       } from '../../../modules/local/splicecounts/mergenonsplitreads'
include { SPLICECOUNTS_COLLECT                  } from '../../../modules/local/splicecounts/collect'
include { FRASER_PSIVALUECALCULATION            } from '../../../modules/local/fraser/psivaluecalculation'
include { FRASER_FILTEREXPRESSION               } from '../../../modules/local/fraser/filterexpression'
include { FRASER_FITHYPERPARAMS                 } from '../../../modules/local/fraser/fithyperparams'
include { FRASER_FITAUTOENCODER                 } from '../../../modules/local/fraser/fitautoencoder'
include { FRASER_ANNOTATEGENES                  } from '../../../modules/local/fraser/annotategenes'
include { FRASER_CALCULATESTATS                 } from '../../../modules/local/fraser/calculatestats'
include { FRASER_EXTRACTRESULTS                 } from '../../../modules/local/fraser/extractresults'

workflow ABERRANTSPLICING {
    take:
    inputs                      // queue channel: [ val(meta), path(bam), path(bai), path(splice_count) ]
    txdb                        // queue channel: The TxDb object for the genome for each gene annotation
    gene_name_mapping           // queue channel: The gene name mapping file in TSV format for each gene annotation
    samplesheet                 // file:          The pipeline samplesheet in TSV format
    fraser_version              // string:        Fraser version to use for aberrant splicing analysis
    include_groups              // list:          A list of groups to include in the aberrant splicing analysis
    genes_to_test               // value channel: [ val(meta), path(genes_to_test) ]
    hpo_file                    // value channel: [ val(meta), path(hpo_file) ]
    aberrant_splicing_config_R  // file:          Path to the R script with configuration for aberrant splicing analysis

    main:
    def ch_versions = Channel.empty()

    def inputs_to_analyse_unfiltered = Channel.empty()
    if (include_groups) {
        // Filter out the BAM files that don't have a group in the include_groups list
        inputs_to_analyse_unfiltered = inputs.filter { it ->
            it[0].drop_group.tokenize(",").intersect(include_groups).size() > 0
        }
    } else {
        inputs_to_analyse_unfiltered = inputs
    }

    // Filter out samples that have external splice counts
    def inputs_to_analyse = inputs_to_analyse_unfiltered
        .filter { _meta, _bam, _bai, splice_counts -> !splice_counts }
        .map { meta, bam, bai, _splice_counts ->
            [ meta, bam, bai ]
        }

    //
    // Create a datasets file
    //

    def ch_group_datasets = inputs_to_analyse
        .map { meta, bam, bai ->
            def groups = include_groups ? meta.drop_group.tokenize(",").intersect(include_groups) : meta.drop_group.tokenize(",")
            [ groups, meta, bam, bai ]
        }
        .transpose(by:0)
        .groupTuple()
        .tap { ch_group_files }
        .map { group, metas, bams, _bais ->
            def header = ["sampleID\tbamFile\tpairedEnd"]
            def lines = []
            def index = 0
            metas.each {
                lines += [
                    metas[index].id,
                    bams[index].name,
                    metas[index].paired_end ? "TRUE" : "FALSE"
                ].join("\t")
                index += 1
            }
            return [group, header + lines.sort()]
        }
        .collectFile(newLine:true, storeDir:"${workflow.workDir}/aberrant_splicing_datasets/") { group, lines ->
            [ "${group}.tsv", lines.join("\n") ]
        }

    // Copy the datasets files to the output directory
    ch_group_datasets.subscribe { dataset ->
            dataset.copyTo("${params.outdir}/processed_data/aberrant_splicing/annotations/${dataset.name}")
        }

    //
    // SPLICECOUNTS_INIT
    //

    def ch_abberant_splicing_input = ch_group_datasets
        .map { dataset ->
            def id = dataset.name.replace(".tsv", "")
            [ id, dataset ]
        }
        .join(ch_group_files, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { id, tsv, metas, bams, bais ->
            def meta = [
                id: id,
                drop_group: id,
                samples: metas.collect { entry -> entry.id }.sort().join(","), // Sort the sample IDs for reproducibility
                group_size: metas.size()
            ]
            dataset: [ meta, tsv ]
            bams: [ meta, bams, bais ]
        }

    SPLICECOUNTS_INIT(
        ch_abberant_splicing_input.dataset.map { meta, dataset ->
            [ meta, dataset, meta.drop_group ]
        },
        fraser_version,
        aberrant_splicing_config_R,
    )
    ch_versions = ch_versions.mix(SPLICECOUNTS_INIT.out.versions.first())

    //
    // SPLICECOUNTS_COUNTSPLITREADS
    //

    def ch_bams_per_sample = inputs_to_analyse
        .map { meta, bam, bai ->
            [ meta.id, bam, bai ]
        }

    // Note for later: Rethink this so that split counting only happens once per sample
    def splitreadssamplewise_input = SPLICECOUNTS_INIT.out.fdsobj
        .map { meta, fds ->
            [ meta.samples.tokenize(","), meta, fds ]
        }
        .transpose(by:0)
        .combine(ch_bams_per_sample, by:0)
        .map { sample, meta, fds, bam, bai ->
            def new_meta = meta + [id:sample]
            [ new_meta, fds, bam, bai, meta.drop_group, sample ]
        }

    SPLICECOUNTS_COUNTSPLITREADS(
        splitreadssamplewise_input,
        params.as_keep_non_standard_chrs,
        params.as_recount,
        params.genome,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(SPLICECOUNTS_COUNTSPLITREADS.out.versions.first())

    //
    // SPLICECOUNTS_MERGESPLITREADS
    //

    def ch_splitreadsmerge_input = SPLICECOUNTS_COUNTSPLITREADS.out.split_counts
        .map { meta, split_counts ->
            def new_meta = meta + [id:meta.drop_group]
            [ groupKey(new_meta, meta.group_size), split_counts ]
        }
        .groupTuple()
        .join(SPLICECOUNTS_INIT.out.fdsobj, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_abberant_splicing_input.bams, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, split_counts, fds, bams, bais ->
            [ meta, fds, split_counts, bams, bais, meta.drop_group ]
        }

    SPLICECOUNTS_MERGESPLITREADS(
        ch_splitreadsmerge_input,
        params.as_min_expression_in_one_sample,
        params.as_recount,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(SPLICECOUNTS_MERGESPLITREADS.out.versions.first())

    //
    // SPLICECOUNTS_COUNTNONSPLITREADS
    //

    def ch_nonsplitreadssamplewise_input = SPLICECOUNTS_MERGESPLITREADS.out.fdsobj
        .join(SPLICECOUNTS_MERGESPLITREADS.out.cache, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, fds, cache ->
            [ meta.samples.tokenize(","), meta, fds, cache ]
        }
        .transpose(by:0)
        .combine(ch_bams_per_sample, by:0)
        .map { sample, meta, fds, cache, bam, bai ->
            def new_meta = meta + [id:sample]
            [ new_meta, fds, cache.resolve("spliceSites_splitCounts.rds"), bam, bai, meta.drop_group, sample ]
        }

    SPLICECOUNTS_COUNTNONSPLITREADS(
        ch_nonsplitreadssamplewise_input,
        params.as_long_read,
        params.as_recount,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(SPLICECOUNTS_COUNTNONSPLITREADS.out.versions.first())

    //
    // SPLICECOUNTS_MERGENONSPLITREADS
    //

    def ch_nonsplitreadsmerge_input = SPLICECOUNTS_COUNTNONSPLITREADS.out.non_split_counts
        .map { meta, non_split_counts ->
            def new_meta = meta + [id:meta.drop_group]
            [ groupKey(new_meta, meta.group_size), non_split_counts ]
        }
        .groupTuple()
        .join(SPLICECOUNTS_MERGESPLITREADS.out.fdsobj, failOnMismatch: true, failOnDuplicate: true)
        .join(SPLICECOUNTS_MERGESPLITREADS.out.cache, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_abberant_splicing_input.bams, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, non_split_counts, fds, cache, bams, bais ->

            [ meta, fds, cache.resolve("gRanges_NonSplitCounts.rds"), non_split_counts, bams, bais, meta.drop_group ]
        }

    SPLICECOUNTS_MERGENONSPLITREADS(
        ch_nonsplitreadsmerge_input,
        params.as_long_read,
        params.as_recount,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(SPLICECOUNTS_MERGENONSPLITREADS.out.versions.first())

    //
    // SPLICECOUNTS_COLLECT
    //

    def ch_collect_input = SPLICECOUNTS_MERGENONSPLITREADS.out.fdsobj
        .join(SPLICECOUNTS_MERGESPLITREADS.out.cache, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, fds, cache ->
            [
                meta,
                fds,
                cache.resolve("gRanges_splitCounts.rds"),
                cache.resolve("spliceSites_splitCounts.rds"),
                meta.drop_group
            ]
        }

    SPLICECOUNTS_COLLECT(
        ch_collect_input,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(SPLICECOUNTS_COLLECT.out.versions.first())

    //
    // FRASER_PSIVALUECALCULATION
    //

    def ch_psivaluecalculation_input = SPLICECOUNTS_COLLECT.out.fdsobj
        .map { meta, fds -> [ meta, fds, meta.drop_group ] }

    FRASER_PSIVALUECALCULATION(
        ch_psivaluecalculation_input,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(FRASER_PSIVALUECALCULATION.out.versions.first())

    //
    // FRASER_FILTEREXPRESSION
    //

    def splice_counts = [:]

    samplesheet.splitCsv(header:true, sep:"\t")
        .each { row ->
            if (!row.SPLICE_COUNTS_DIR) {
                return
            }
            def groups = row.DROP_GROUP.tokenize(",")
            groups.each { group ->
                splice_counts[group] = splice_counts.get(group, []) + [[id:row.RNA_ID, file:row.SPLICE_COUNTS_DIR]]
            }
        }

    def ch_filterexpression_input = FRASER_PSIVALUECALCULATION.out.fdsobj
        .map { meta, fds ->
            def splice_counts_group = splice_counts.get(meta.drop_group, [])
            def count_dirs = splice_counts_group.collect { it.file }.unique() // Prevent the same directory from being used multiple times
            def count_ids = splice_counts_group.collect { it.id }
            def new_meta = meta + [group_size: meta.group_size + count_ids.size(), samples: (meta.samples.tokenize(",") + count_ids).sort().join(",")]
            [ new_meta, fds, count_dirs, count_ids, meta.drop_group ]
        }

    FRASER_FILTEREXPRESSION(
        ch_filterexpression_input,
        Channel.value([[id:'samplesheet'], samplesheet]),
        params.as_min_expression_in_one_sample,
        params.as_quantile_for_filtering,
        params.as_quantile_min_expression,
        params.as_min_delta_psi,
        !params.as_skip_filter,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(FRASER_FILTEREXPRESSION.out.versions.first())

    //
    // FRASER_FITHYPERPARAMS
    //

    def ch_fithyperparameters_input = FRASER_FILTEREXPRESSION.out.fdsobj
        .map { meta, fds -> [ meta, fds, meta.drop_group ] }

    FRASER_FITHYPERPARAMS(
        ch_fithyperparameters_input,
        params.random_seed,
        params.as_implementation,
        params.as_max_tested_dimension_proportion,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(FRASER_FITHYPERPARAMS.out.versions.first())

    //
    // FRASER_FITAUTOENCODER
    //

    def ch_fitautoencoder_input = FRASER_FITHYPERPARAMS.out.fdsobj
        .map { meta, fds -> [ meta, fds, meta.drop_group ] }

    FRASER_FITAUTOENCODER(
        ch_fitautoencoder_input,
        params.as_implementation,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(FRASER_FITAUTOENCODER.out.versions.first())

    //
    // FRASER_ANNOTATEGENES
    //

    def ch_gene_annotations = txdb.join(gene_name_mapping, failOnMismatch: true, failOnDuplicate: true)

    def ch_annotategenes_input = FRASER_FITAUTOENCODER.out.fdsobj
        .combine(ch_gene_annotations)
        .map { meta, fds, meta_annotations, txdb_, gene_name_mapping_ ->
            def new_meta = meta + [id:"${meta.id}.${meta_annotations.id}", annotation:meta_annotations.id]
            [ new_meta, fds, txdb_, gene_name_mapping_, meta.drop_group, meta_annotations.id ]
        }

    FRASER_ANNOTATEGENES(
        ch_annotategenes_input,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(FRASER_ANNOTATEGENES.out.versions.first())

    //
    // FRASER_CALCULATESTATS
    //

    def ch_calculationstats_input = FRASER_ANNOTATEGENES.out.fdsobj
        .map { meta, fds ->
            [ meta, fds, meta.drop_group, meta.annotation, meta.samples.tokenize(",") ]
        }

    FRASER_CALCULATESTATS(
        ch_calculationstats_input,
        genes_to_test,
        fraser_version,
        aberrant_splicing_config_R,
        file("${projectDir}/assets/helpers/parse_subsets_for_FDR.R", checkIfExists: true),
    )
    ch_versions = ch_versions.mix(FRASER_CALCULATESTATS.out.versions.first())

    //
    // FRASER_EXTRACTRESULTS
    //

    def ch_extractresults_input = FRASER_CALCULATESTATS.out.fdsobj
        .map { meta, fds ->
            [meta.annotation, meta, fds]
        }
        .combine(txdb.map { meta, txdb_ -> [meta.id, txdb_] }, by:0)
        .map { annotation, meta, fds, txdb_ ->
            [meta, fds, txdb_, meta.drop_group, annotation, meta.samples.tokenize(",")]
        }

    FRASER_EXTRACTRESULTS(
        ch_extractresults_input,
        Channel.value([[id:'samplesheet'], samplesheet]),
        hpo_file,
        params.as_padj_cutoff,
        params.as_delta_psi_cutoff,
        params.genome,
        fraser_version,
        aberrant_splicing_config_R,
        file("${projectDir}/assets/helpers/add_HPO_cols.R", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(FRASER_EXTRACTRESULTS.out.versions.first())

    emit:
    versions = ch_versions
}
