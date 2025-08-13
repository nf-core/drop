#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/drop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/drop
    Website: https://nf-co.re/drop
    Slack  : https://nfcore.slack.com/channels/drop
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DROP                              } from './workflows/drop'

include { GATK4_CREATESEQUENCEDICTIONARY as GATK4_CREATESEQUENCEDICTIONARY_UCSC } from './modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK4_CREATESEQUENCEDICTIONARY_NCBI } from './modules/nf-core/gatk4/createsequencedictionary/main'

include { PIPELINE_INITIALISATION           } from './subworkflows/local/utils_nfcore_drop_pipeline'
include { PIPELINE_COMPLETION               } from './subworkflows/local/utils_nfcore_drop_pipeline'
include { getGenomeAttribute                } from './subworkflows/local/utils_nfcore_drop_pipeline'

include { samplesheetToList                 } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
// params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // Parse and convert parameters to their expected input format
    //

    if(!params.ucsc_fasta && !params.ncbi_fasta) {
        error "You must provide either a UCSC FASTA file (--ucsc_fasta) or a NCBI FASTA file (--ncbi_fasta)."
    }

    // Use the UCSC FASTA and FAI files if provided, otherwise use the NCBI FASTA and FAI files
    def fasta = params.ucsc_fasta ? Channel.value([[id: 'fasta'], file(params.ucsc_fasta)]) :
        params.ncbi_fasta ? Channel.value([[id: 'fasta'], file(params.ncbi_fasta)]) : Channel.empty()
    def fai = params.ucsc_fai ? Channel.value([[id: 'fai'], file(params.ucsc_fai)]) :
        params.ncbi_fai ? Channel.value([[id: 'fai'], file(params.ncbi_fai)]) : Channel.empty()

    def ucsc_fasta = params.ucsc_fasta ? Channel.value([[id: 'ucsc'], file(params.ucsc_fasta)]) : Channel.of([[id:'ucsc'], []])
    def ucsc_fai = params.ucsc_fai ? Channel.value([[id: 'ucsc'], file(params.ucsc_fai)]) : Channel.of([[id:'ucsc'], []])
    def ncbi_fasta = params.ncbi_fasta ? Channel.value([[id: 'ncbi'], file(params.ncbi_fasta)]) : Channel.of([[id:'ncbi'], []])
    def ncbi_fai = params.ncbi_fai ? Channel.value([[id: 'ncbi'], file(params.ncbi_fai)]) : Channel.of([[id:'ncbi'], []])

    def qc_vcf = params.mae_qc_vcf ?
        Channel.value([[id: 'qc_vcf'], file(params.mae_qc_vcf), params.mae_qc_vcf_tbi ? file(params.mae_qc_vcf_tbi) : []]) :
        [[:], [], []]

    def ucsc_dict = Channel.of([[id:'ucsc'], []])
    if (params.ucsc_dict) {
        ucsc_dict = Channel.value([[id: 'ucsc'], file(params.ucsc_dict)])
    } else if (params.ucsc_fasta) {
        GATK4_CREATESEQUENCEDICTIONARY_UCSC(
            ucsc_fasta
        )
        ucsc_dict = GATK4_CREATESEQUENCEDICTIONARY_UCSC.out.dict.collect()
    }

    def ncbi_dict = Channel.of([[id:'ncbi'], []])
    if (params.ncbi_dict) {
        ncbi_dict = Channel.value([[id: 'ncbi'], file(params.ncbi_dict)])
    } else if (params.ncbi_fasta) {
        GATK4_CREATESEQUENCEDICTIONARY_NCBI(
            ncbi_fasta
        )
        ncbi_dict = GATK4_CREATESEQUENCEDICTIONARY_NCBI.out.dict.collect()
    }

    def samplesheet_file = Channel.value([[id: 'samplesheet'], file(params.input)])

    def hpo_file = params.hpo_file ? Channel.value([[id: 'hpo'], file(params.hpo_file)]) : [[:], []]

    def ec_gene_annotations = params.ec_gene_annotations ? params.ec_gene_annotations.tokenize(",") : []
    def ec_exclude_groups = params.ec_exclude_groups ? params.ec_exclude_groups.tokenize(",") : []

    def ae_groups = params.ae_groups ? params.ae_groups.tokenize(",") : []
    def ae_genes_to_test = params.ae_genes_to_test ? Channel.value([[id: 'genes_to_test'], file(params.ae_genes_to_test)]) : [[:], []]

    def as_groups = params.as_groups ? params.as_groups.tokenize(",") : []
    def as_genes_to_test = params.as_genes_to_test ? Channel.value([[id: 'genes_to_test'], file(params.as_genes_to_test)]) : [[:], []]

    def mae_groups = params.mae_groups ? params.mae_groups.tokenize(",") : []
    def mae_qc_groups = params.mae_qc_groups ? params.mae_qc_groups.tokenize(",") : []
    //
    // WORKFLOW: Run main workflow
    //
    DROP (
        // Global parameters
        PIPELINE_INITIALISATION.out.samplesheet, // derived from --input
        samplesheet_file,
        params.project_title,
        fasta,
        fai,
        ucsc_fasta,
        ucsc_fai,
        ncbi_fasta,
        ncbi_fai,
        ucsc_dict,
        ncbi_dict,
        PIPELINE_INITIALISATION.out.gene_annotation,
        hpo_file,
        qc_vcf,

        // Export counts parameters
        ec_gene_annotations,
        ec_exclude_groups,

        // Aberrant expression parameters
        params.ae_run,
        ae_groups,
        ae_genes_to_test,

        // Aberrant splicing parameters
        params.as_run,
        as_groups,
        params.as_fraser_version,
        as_genes_to_test,

        // Mono Allelic Expression parameters
        params.mae_run,
        mae_groups,
        mae_qc_groups
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        DROP.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
