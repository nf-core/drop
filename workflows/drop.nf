/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_INDEX            } from '../modules/nf-core/samtools/index/main'
include { TABIX_TABIX               } from '../modules/nf-core/tabix/tabix/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'

include { PREPROCESSGENEANNOTATION  } from '../modules/local/preprocessgeneannotation/main'

include { ABERRANTEXPRESSION        } from '../subworkflows/local/aberrantexpression/main'

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_drop_pipeline'

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
    project_title       // string:        title of the project to add to the HTML file
    fasta               // value channel: [ val(meta), path(fasta) ]
    fai                 // value channel: [ val(meta), path(fai) ]
    gene_annotation     // queue channel: [ val(meta), path(gtf) ]
    hpo_file            // value channel: [ val(meta), path(hpo_file) ]

    // Export count parameters
    ec_gene_annotations // list:          A list of gene annotations to export the counts of
    ec_exclude_groups   // list:          A list of groups to exclude from the counts export

    // Aberrant expression parameters
    ae_run              // boolean:       Run aberrant expression analysis
    ae_groups           // list:          A list of groups to exclude from the aberrant expression analysis
    ae_genes_to_test    // map:           A map containing the names of genes to test

    // Aberrant splicing parameters
    as_run              // boolean:       Run aberrant splicing analysis
    as_groups           // list:          A list of groups to exclude from the aberrant splicing analysis
    as_fraser_version   // string:        Fraser version to use for aberrant splicing analysis
    as_genes_to_test    // map:           A map containing the names of genes to test

    // Mono Allelic Expression parameters
    mae_run             // boolean:       Run mono allelic expression analysis
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

    def bams_to_index = preprocess.bam.branch { meta, bam, bai ->
        yes: !bai && bam
            return [ meta, bam ]
        no: true
    }

    SAMTOOLS_INDEX(
        bams_to_index.yes
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def bams = bams_to_index.yes
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .mix(bams_to_index.no)

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
            // TODO: Create channels for each subworkflow here
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

    if(ae_run) {
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
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'drop_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
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
