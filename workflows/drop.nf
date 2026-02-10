/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_INDEX                            } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CONVERT_UCSC } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CONVERT_NCBI } from '../modules/nf-core/samtools/convert/main'
include { TABIX_TABIX                               } from '../modules/nf-core/tabix/tabix/main'
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'

include { PREPROCESSGENEANNOTATION                  } from '../modules/local/preprocessgeneannotation/main'
include { ABERRANTEXPRESSION                        } from '../subworkflows/local/aberrantexpression/main'
include { ABERRANTSPLICING                          } from '../subworkflows/local/aberrantsplicing/main'
include { MAE                                       } from '../subworkflows/local/mae/main'

include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_drop_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DROP {

    take:
    // Global parameters
    samplesheet         // queue channel: samplesheet read in from --input
    samplesheet_file    // value channel: [ val(meta), path(samplesheet) ]
    ucsc_fasta          // value channel: [ val(meta), path(ucsc_fasta) ]
    ucsc_fai            // value channel: [ val(meta), path(ucsc_fai) ]
    ncbi_fasta          // value channel: [ val(meta), path(ncbi_fasta) ]
    ncbi_fai            // value channel: [ val(meta), path(ncbi_fai) ]
    ucsc_dict           // value channel: [ val(meta), path(ucsc_dict) ]
    ncbi_dict           // value channel: [ val(meta), path(ncbi_dict) ]
    gene_annotation     // queue channel: [ val(meta), path(gtf) ]
    hpo_file            // value channel: [ val(meta), path(hpo_file) ]
    qc_vcf              // value channel: [ val(meta), path(vcf), path(tbi) ]

    // Aberrant expression parameters
    ae_skip             // boolean:       Skip aberrant expression analysis
    ae_groups           // list:          A list of groups to exclude from the aberrant expression analysis
    ae_genes_to_test    // map:           A map containing the names of genes to test

    // Aberrant splicing parameters
    as_skip             // boolean:       Skip aberrant splicing analysis
    as_groups           // list:          A list of groups to exclude from the aberrant splicing analysis
    as_fraser_version   // string:        Fraser version to use for aberrant splicing analysis
    as_genes_to_test    // map:           A map containing the names of genes to test

    // Mono Allelic Expression parameters
    mae_skip            // boolean:       Skip mono allelic expression analysis
    mae_groups          // list:          A list of groups to exclude from the mono allelic expression analysis
    mae_qc_groups       // list:          A list of groups to exclude from QC steps in the mono allelic expression analysis

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    def preprocess = samplesheet.multiMap { meta, rna_bam, rna_bai, dna_vcf, dna_tbi, gene_counts, splice_counts ->
        bam: [ meta, rna_bam, rna_bai ]
        vcf: [ meta, dna_vcf, dna_tbi ]
        counts: [ meta, gene_counts, splice_counts ]
    }

    //
    // Preprocess BAM files
    //

    def bams_branch = preprocess.bam.branch { meta, bam, bai ->
        cram_ucsc: bam.extension == "cram" && meta.genome == "ucsc"
            return [ meta, bam, bai ]
        cram_ncbi: bam.extension == "cram" && meta.genome == "ncbi"
            return [ meta, bam, bai ]
        to_index: !bai && bam
            return [ meta, bam ]
        indexed: true
    }

    SAMTOOLS_CONVERT_UCSC(
        bams_branch.cram_ucsc,
        ucsc_fasta,
        ucsc_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT_UCSC.out.versions.first())
    def ch_ucsc_converted_bams = SAMTOOLS_CONVERT_UCSC.out.bam.join(SAMTOOLS_CONVERT_UCSC.out.bai, failOnDuplicate:true, failOnMismatch:true)

    SAMTOOLS_CONVERT_NCBI(
        bams_branch.cram_ncbi,
        ncbi_fasta,
        ncbi_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT_NCBI.out.versions.first())
    def ch_ncbi_converted_bams = SAMTOOLS_CONVERT_NCBI.out.bam.join(SAMTOOLS_CONVERT_NCBI.out.bai, failOnDuplicate:true, failOnMismatch:true)

    SAMTOOLS_INDEX(
        bams_branch.to_index
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def bams = bams_branch.to_index
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_ucsc_converted_bams)
        .mix(ch_ncbi_converted_bams)
        .mix(bams_branch.indexed)

    //
    // Preprocess VCF files
    //

    def vcfs_to_index = preprocess.vcf.branch { meta, vcf, tbi ->
        yes: vcf && !tbi
            return [ meta, vcf ]
        no: true
    }

    TABIX_TABIX(
        vcfs_to_index.yes
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    def vcfs = vcfs_to_index.yes
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .mix(vcfs_to_index.no)

    //
    // Define inputs for each subworkflow
    //

    def input = bams
        .join(vcfs, failOnDuplicate:true, failOnMismatch:true)
        .join(preprocess.counts, failOnDuplicate:true, failOnMismatch:true)
        .multiMap { meta, rna_bam, rna_bai, dna_vcf, dna_tbi, gene_counts, splice_counts ->
            abberantexpression: [ meta, rna_bam, rna_bai, gene_counts ]
            aberrantsplicing: [ meta, rna_bam, rna_bai, splice_counts ]
            mae: [ meta, dna_vcf, dna_tbi, rna_bam, rna_bai ]
        }

    //
    // Preprocess gene annotation
    //

    PREPROCESSGENEANNOTATION(
        gene_annotation
    )
    ch_versions = ch_versions.mix(PREPROCESSGENEANNOTATION.out.versions.first())

    //
    // Abberant expression
    //

    if(!ae_skip) {
        ABERRANTEXPRESSION(
            input.abberantexpression,
            PREPROCESSGENEANNOTATION.out.count_ranges,
            PREPROCESSGENEANNOTATION.out.txdb,
            PREPROCESSGENEANNOTATION.out.gene_name_mapping,
            samplesheet_file,
            ae_genes_to_test,
            hpo_file,
            ae_groups,
        )
        ch_versions = ch_versions.mix(ABERRANTEXPRESSION.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ABERRANTEXPRESSION.out.count_report)
        ch_multiqc_files = ch_multiqc_files.mix(ABERRANTEXPRESSION.out.outrider_report)

    }

    //
    // Aberrant splicing
    //

    if (!as_skip) {
        ABERRANTSPLICING(
            input.aberrantsplicing,
            PREPROCESSGENEANNOTATION.out.txdb,
            PREPROCESSGENEANNOTATION.out.gene_name_mapping,
            file(params.input),
            as_fraser_version,
            as_groups,
            as_genes_to_test,
            hpo_file,
            file("${projectDir}/assets/helpers/aberrant_splicing_config.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(ABERRANTSPLICING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ABERRANTSPLICING.out.count_report)
        ch_multiqc_files = ch_multiqc_files.mix(ABERRANTSPLICING.out.fraser_report)
    }

    //
    // Mono Allelic Expression
    //

    if(!mae_skip) {
        MAE(
            input.mae,
            ucsc_fasta,
            ucsc_fai,
            ucsc_dict,
            ncbi_fasta,
            ncbi_fai,
            ncbi_dict,
            gene_annotation,
            samplesheet_file,
            qc_vcf,
            mae_groups,
            mae_qc_groups,
            Channel.value(file("${projectDir}/assets/chr_NCBI_UCSC.txt", checkIfExists: true)),
            Channel.value(file("${projectDir}/assets/chr_UCSC_NCBI.txt", checkIfExists: true))
        )
        ch_versions = ch_versions.mix(MAE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(MAE.out.mae_report)
        ch_multiqc_files = ch_multiqc_files.mix(MAE.out.maeqc_report)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'drop_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
