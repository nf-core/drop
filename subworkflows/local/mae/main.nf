include { MAE_CREATESNVS                } from '../../../modules/local/mae/createsnvs/main'
include { MAE_ALLELICCOUNTS             } from '../../../modules/local/mae/alleliccounts/main'
include { MAE_DESEQ                     } from '../../../modules/local/mae/deseq/main'
include { MAE_GENENAMEMAPPING           } from '../../../modules/local/mae/genenamemapping/main'
include { MAE_RESULTS                   } from '../../../modules/local/mae/results/main'
include { MAEQC_DNARNADESEQ             } from '../../../modules/local/maeqc/dnarnadeseq/main'
include { MAEQC_DNARNAMATRIX            } from '../../../modules/local/maeqc/dnarnamatrix/main'
include { MAEQC_DNARNAMATRIXPLOT        } from '../../../modules/local/maeqc/dnarnamatrixplot/main'
include { MULTIQC as MULTIQC_MAE        } from '../../../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_MAEQC      } from '../../../modules/nf-core/multiqc/main'
include { STRIP_NAV_MULTIQC as STRIP_NAV_MAE     } from '../../../modules/local/strip_nav_multiqc/'
include { STRIP_NAV_MULTIQC as STRIP_NAV_MAEQC       } from '../../../modules/local/strip_nav_multiqc/'

workflow MAE {
    take:
    input               // queue channel: [ val(meta), path(vcf), path(tbi), path(bam), path(bai) ]
    ucsc_fasta          // value channel: [ val(meta), path(ucsc_fasta) ]
    ucsc_fai            // value channel: [ val(meta), path(ucsc_fai) ]
    ucsc_dict           // value channel: [ val(meta), path(ucsc_dict) ]
    ncbi_fasta          // value channel: [ val(meta), path(ncbi_fasta) ]
    ncbi_fai            // value channel: [ val(meta), path(ncbi_fai) ]
    ncbi_dict           // value channel: [ val(meta), path(ncbi_dict) ]
    gene_annotation     // queue channel: [ val(meta), path(gtf) ]
    samplesheet         // value channel: [ val(meta), path(samplesheet) ]
    qc_vcf              // value channel: [ val(meta), path(vcf), path(tbi) ]
    include_groups      // list         : A list of groups to include in the mono allelic expression analysis
    include_qc_groups   // list         : A list of groups to include in the QC steps of the mono allelic expression analysis
    ncbi2ucsc           // value channel: path to the NCBI to UCSC mapping file
    ucsc2ncbi           // value channel: path to the UCSC to NCBI mapping file

    main:
    def ch_versions = Channel.empty()

    def ch_ucsc_refs = ucsc_fasta
        .join(ucsc_fai, failOnDuplicate:true, failOnMismatch:true)
        .join(ucsc_dict, failOnDuplicate:true, failOnMismatch:true)

    def ch_ncbi_refs = ncbi_fasta
        .join(ncbi_fai, failOnDuplicate:true, failOnMismatch:true)
        .join(ncbi_dict, failOnDuplicate:true, failOnMismatch:true)

    def ch_references = ch_ucsc_refs.mix(ch_ncbi_refs)

    def ch_filtered_inputs = Channel.empty()
    if (include_groups) {
        ch_filtered_inputs = input.filter { meta, _vcf, _tbi, _bam, _bai ->
            meta.drop_group.tokenize(",").intersect(include_groups).size() > 0
        }
        .map { meta, vcf, tbi, bam, bai ->
            def new_meta = meta + [drop_group: meta.drop_group.tokenize(",").intersect(include_groups).join(",")]
            [ new_meta, vcf, tbi, bam, bai ]
        }
    } else {
        ch_filtered_inputs = input
    }

    def ch_createsnvs_input = ch_filtered_inputs
        .filter { _meta, vcf, tbi, bam, bai ->
            vcf && tbi && bam && bai
        }
        .map { meta, vcf, tbi, bam, bai ->
            [ meta, vcf, tbi, bam, bai, meta.dna_id ]
        }

    MAE_CREATESNVS(
        ch_createsnvs_input,
        ncbi2ucsc,
        ucsc2ncbi
    )
    ch_versions = ch_versions.mix(MAE_CREATESNVS.out.versions.first())

    def ch_alleliccounts_input = MAE_CREATESNVS.out.vcf
        .join(MAE_CREATESNVS.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .join(ch_createsnvs_input, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, new_vcf, new_tbi, _old_vcf, _old_tbi, bam, bai, id ->
            [ [id:meta.genome], meta, new_vcf, new_tbi, bam, bai, id ]
        }
        .combine(ch_references, by:0)
        .map { _ref_meta, meta, vcf, tbi, bam, bai, id, fasta, fai, dict ->
            [ meta, vcf, tbi, bam, bai, id, fasta, fai, dict ]
        }

    MAE_ALLELICCOUNTS(
        ch_alleliccounts_input,
        !params.mae_gatk_header_check,
        ncbi2ucsc,
        ucsc2ncbi
    )
    ch_versions = ch_versions.mix(MAE_ALLELICCOUNTS.out.versions.first())

    MAE_DESEQ(
        MAE_ALLELICCOUNTS.out.csv,
        params.mae_add_af,
        params.genome,
        params.mae_max_af
    )
    ch_versions = ch_versions.mix(MAE_DESEQ.out.versions.first())

    MAE_GENENAMEMAPPING(
        gene_annotation
    )
    ch_versions = ch_versions.mix(MAE_GENENAMEMAPPING.out.versions.first())

    def ch_results_input = MAE_DESEQ.out.res
        .combine(MAE_GENENAMEMAPPING.out.tsv)
        .filter { meta, _res, annotation_meta, _tsv ->
            // Filter out the gene annotations that should not be used
            meta.gene_annotation != "" ? meta.gene_annotation == annotation_meta.id : true
        }
        .map { meta, res, annotations_meta, tsv ->
            def new_meta = [
                id: meta.drop_group,
                drop_group: meta.drop_group,
                annotation: annotations_meta.id
            ]
            [ groupKey(new_meta, meta.drop_group_ann_counts.get(meta.drop_group).get(annotations_meta.id)), res, tsv ]
        }
        .groupTuple()
        .map { meta, res_files, tsv_files ->
            // Make sure the gene name mapping TSV is only provided once
            [ meta, res_files, tsv_files.find { it -> it } ]
        }

    MAE_RESULTS(
        ch_results_input,
        params.mae_allelic_ratio_cutoff,
        params.mae_padj_cutoff,
        params.mae_max_var_freq_cohort
    )
    ch_versions = ch_versions.mix(MAE_RESULTS.out.versions.first())

    // TODO implement mae_qc_groups logic
    MAEQC_DNARNADESEQ(
        MAE_ALLELICCOUNTS.out.csv
    )
    ch_versions = ch_versions.mix(MAEQC_DNARNADESEQ.out.versions.first())

    def ch_creatematrixdnarnacor_input = MAEQC_DNARNADESEQ.out.ods_with_pvals
        .join(MAE_CREATESNVS.out.vcf, failOnDuplicate:true, failOnMismatch:true)
        .join(MAE_CREATESNVS.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, ods, vcf, tbi ->
            [ meta.drop_group.tokenize(","), meta, ods, vcf, tbi ]
        }
        .transpose(by:0)
        .map { group, meta, ods, vcf, tbi ->
            def new_meta = [id:group, drop_group:group]
            [ groupKey(new_meta, meta.drop_group_counts.get(group)), ods, vcf, tbi ]
        }
        .groupTuple()
        .combine(samplesheet)
        .map { meta, ods, vcfs, tbis, _ss_meta, samplesheet_file ->
            // Sort VCFs and IDs to make sure the order matches
            [ meta, ods, vcfs.sort { file -> file.name }, tbis, samplesheet_file, meta.drop_group ]
        }

    MAEQC_DNARNAMATRIX(
        ch_creatematrixdnarnacor_input,
        qc_vcf
    )
    ch_versions = ch_versions.mix(MAEQC_DNARNAMATRIX.out.versions.first())

    ch_dnarnamatrixplot_input = MAEQC_DNARNAMATRIX.out.mat_qc
        .combine(samplesheet)
        .map { meta, mat_qc, _ss_meta, samplesheet_file ->
            [ meta, mat_qc, samplesheet_file, meta.drop_group ]
        }
        .view {println "[DEBUG ch_dnarnamatrixplot_input] $it"}

    MAEQC_DNARNAMATRIXPLOT(
        ch_dnarnamatrixplot_input,
        params.mae_dna_rna_match_cutoff
    )
    ch_versions = ch_versions.mix(MAEQC_DNARNAMATRIXPLOT.out.versions.first())

    //
    // multiqc for mae
    //
    def multiqc_mae_input = MAE_RESULTS.out.mae_overview
        .join(MAE_RESULTS.out.cascade_plot)
        .join(MAE_RESULTS.out.variant_frequency)
        .join(MAE_RESULTS.out.median_of_each_category)
        .join(MAE_RESULTS.out.results)
        .map { tup ->
        def meta  = tup[0]
        def files = tup.drop(1).findAll { it instanceof Path } // drop meta and null
        def tag   = "${meta.id}" // tag for publishDir
        [ tag, files ]
    }
    .tap { mae_by_tag }

    def ch_extra_cfg_mae = mae_by_tag
        .collectFile { tag, _ ->
            def yaml = """
            output_fn_name: "[TAG:mae_${tag}]_multiqc_report.html"
            data_dir_name:  "[TAG:mae_${tag}]_multiqc_data"
            plots_dir_name: "[TAG:mae_${tag}]_multiqc_plots"
            """.stripIndent()
            [ "[TAG:mae_${tag}]_config.yml", yaml ]
        }
        .map { Path cfg ->
        def m = (cfg.getName() =~ /\[TAG:mae_([^\]]+)\]_config\.yml$/)
        def tag = m ? m[0][1] : cfg.baseName
        [ tag, cfg ]
    }

    def ch_mae_bundle = multiqc_mae_input
        .join(ch_extra_cfg_mae)
        .map { _tag, files, cfg -> [files, cfg] }

    MULTIQC_MAE(
        ch_mae_bundle.map { it[0] },
        Channel.fromPath("$projectDir/assets/multiqc_configs/multiqc_mae_config.yml", checkIfExists: true).collect(),
        ch_mae_bundle.map { it[1] },
        [],
        [],
        []
    )

    //
    // multiqc for maeqc
    //
    def multiqc_maeqc_input = MAEQC_DNARNAMATRIXPLOT.out.matching_values_distribution
        .join(MAEQC_DNARNAMATRIXPLOT.out.heatmap_matching_variants)
        .join(MAEQC_DNARNAMATRIXPLOT.out.identify_matching_samples)
        .join(MAEQC_DNARNAMATRIXPLOT.out.false_matches)
        .join(MAEQC_DNARNAMATRIXPLOT.out.false_mismatches)
        .map { tup ->
        def meta  = tup[0]
        def files = tup.drop(1).findAll { it instanceof Path } // drop meta and null
        def tag   = "${meta.id}" // tag for publishDir
        [ tag, files ]
    }
    .tap { maeqc_by_tag }

    def ch_extra_cfg_maeqc = maeqc_by_tag
        .collectFile { tag, _ ->
            def yaml = """
            output_fn_name: "[TAG:maeqc_${tag}]_multiqc_report.html"
            data_dir_name:  "[TAG:maeqc_${tag}]_multiqc_data"
            plots_dir_name: "[TAG:maeqc_${tag}]_multiqc_plots"
            """.stripIndent()
            [ "[TAG:maeqc_${tag}]_config.yml", yaml ]
        }
        .map { Path cfg ->
        def m = (cfg.getName() =~ /\[TAG:maeqc_([^\]]+)\]_config\.yml$/)
        def tag = m ? m[0][1] : cfg.baseName
        [ tag, cfg ]
    }

    def ch_maeqc_bundle = multiqc_maeqc_input
        .join(ch_extra_cfg_maeqc)
        .map { _tag, files, cfg -> [files, cfg] }

    MULTIQC_MAEQC(
        ch_maeqc_bundle.map { it[0] },
        Channel.fromPath("$projectDir/assets/multiqc_configs/multiqc_maeqc_config.yml", checkIfExists: true).collect(),
        ch_maeqc_bundle.map { it[1] },
        [],
        [],
        []
    )

    STRIP_NAV_MAE(MULTIQC_MAE.out.report)
    STRIP_NAV_MAEQC(MULTIQC_MAEQC.out.report)

    emit:
    versions = ch_versions
    mae_report = STRIP_NAV_MAE.out.report.toList()
    maeqc_report = STRIP_NAV_MAEQC.out.report.toList()
}
