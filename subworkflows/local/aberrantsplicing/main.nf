include { COUNTRNA_INIT                     } from '../../../modules/local/countrna/init'
include { COUNTRNA_SPLITREADSSAMPLEWISE     } from '../../../modules/local/countrna/splitreadssamplewise'
include { COUNTRNA_SPLITREADSMERGE          } from '../../../modules/local/countrna/splitreadsmerge'
include { COUNTRNA_NONSPLITREADSSAMPLEWISE  } from '../../../modules/local/countrna/nonsplitreadssamplewise'
include { COUNTRNA_NONSPLITREADSMERGE       } from '../../../modules/local/countrna/nonsplitreadsmerge'

workflow ABERRANTSPLICING {
    take:
    inputs                      // queue channel: [ val(meta), path(bam), path(bai), path(vcf), path(tbi), path(gene_count), path(splice_count) ]
    fraser_version              // string:        Fraser version to use for aberrant splicing analysis
    include_groups              // list:          A list of groups to include in the aberrant splicing analysis
    genes_to_test               // value channel: [ val(meta), path(genes_to_test) ]
    aberrant_splicing_config_R  // file:          Path to the R script with configuration for aberrant splicing analysis

    main:
    def ch_versions = Channel.empty()

    def inputs_to_analyse = Channel.empty()
    if (include_groups) {
        // Filter out the BAM files that don't have a group in the include_groups list
        inputs_to_analyse = inputs.filter { it ->
            it[0].drop_group.tokenize(",").intersect(include_groups).size() > 0
        }
    } else {
        inputs_to_analyse = inputs
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
        .collectFile(newLine:true, storeDir:"${params.outdir}/processed_data/aberrant_splicing/annotations/") { group, lines ->
            [ "${group}.tsv", lines.join("\n") ]
        }

    //
    // COUNTRNA_INIT
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

    COUNTRNA_INIT(
        ch_abberant_splicing_input.dataset.map { meta, dataset ->
            [ meta, dataset, meta.drop_group ]
        },
        fraser_version,
        aberrant_splicing_config_R,
    )
    ch_versions = ch_versions.mix(COUNTRNA_INIT.out.versions.first())

    //
    // COUNTRNA_SPLITREADSSAMPLEWISE
    //

    def ch_bams_per_sample = inputs_to_analyse
        .map { meta, bam, bai ->
            [ meta.id, bam, bai ]
        }

    // Note for later: Rethink this so that split counting only happens once per sample
    def splitreadssamplewise_input = COUNTRNA_INIT.out.fdsobj
        .map { meta, fds ->
            [ meta.samples.tokenize(","), meta, fds ]
        }
        .transpose(by:0)
        .combine(ch_bams_per_sample, by:0)
        .map { sample, meta, fds, bam, bai ->
            def new_meta = meta + [id:sample]
            [ new_meta, fds, bam, bai, meta.drop_group, sample ]
        }

    COUNTRNA_SPLITREADSSAMPLEWISE(
        splitreadssamplewise_input,
        params.as_keep_non_standard_chrs,
        params.as_recount,
        params.genome,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(COUNTRNA_SPLITREADSSAMPLEWISE.out.versions.first())

    //
    // COUNTRNA_SPLITREADSMERGE
    //

    def ch_splitreadsmerge_input = COUNTRNA_SPLITREADSSAMPLEWISE.out.split_counts
        .map { meta, split_counts ->
            def new_meta = meta + [id:meta.drop_group]
            [ groupKey(new_meta, meta.group_size), split_counts ]
        }
        .groupTuple()
        .join(COUNTRNA_INIT.out.fdsobj, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_abberant_splicing_input.bams, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, split_counts, fds, bams, bais ->
            [ meta, fds, split_counts, bams, bais, meta.drop_group ]
        }

    COUNTRNA_SPLITREADSMERGE(
        ch_splitreadsmerge_input,
        params.as_min_expression_in_one_sample,
        params.as_recount,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(COUNTRNA_SPLITREADSMERGE.out.versions.first())

    //
    // COUNTRNA_NONSPLITREADSSAMPLEWISE
    //

    def ch_nonsplitreadssamplewise_input = COUNTRNA_SPLITREADSMERGE.out.fdsobj
        .join(COUNTRNA_SPLITREADSMERGE.out.cache, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, fds, cache ->
            [ meta.samples.tokenize(","), meta, fds, cache ]
        }
        .transpose(by:0)
        .combine(ch_bams_per_sample, by:0)
        .map { sample, meta, fds, cache, bam, bai ->
            def new_meta = meta + [id:sample]
            [ new_meta, fds, cache.resolve("raw-local-${meta.drop_group}/spliceSites_splitCounts.rds"), bam, bai, meta.drop_group, sample ]
        }

    COUNTRNA_NONSPLITREADSSAMPLEWISE(
        ch_nonsplitreadssamplewise_input,
        params.as_long_read,
        params.as_recount,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(COUNTRNA_NONSPLITREADSSAMPLEWISE.out.versions.first())

    //
    // COUNTRNA_NONSPLITREADSMERGE
    //

    def ch_nonsplitreadsmerge_input = COUNTRNA_NONSPLITREADSSAMPLEWISE.out.non_split_counts
        .map { meta, non_split_counts ->
            def new_meta = meta + [id:meta.drop_group]
            [ groupKey(new_meta, meta.group_size), non_split_counts ]
        }
        .groupTuple()
        .join(COUNTRNA_SPLITREADSMERGE.out.fdsobj, failOnMismatch: true, failOnDuplicate: true)
        .join(COUNTRNA_SPLITREADSMERGE.out.cache, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_abberant_splicing_input.bams, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, non_split_counts, fds, cache, bams, bais ->

            [ meta, fds, cache.resolve("raw-local-${meta.drop_group}/gRanges_NonSplitCounts.rds"), non_split_counts, bams, bais, meta.drop_group ]
        }

    COUNTRNA_NONSPLITREADSMERGE(
        ch_nonsplitreadsmerge_input,
        params.as_long_read,
        params.as_recount,
        fraser_version,
        aberrant_splicing_config_R
    )
    ch_versions = ch_versions.mix(COUNTRNA_NONSPLITREADSMERGE.out.versions.first())

    COUNTRNA_NONSPLITREADSMERGE.out.cache.view()

    emit:
    versions = ch_versions
}
