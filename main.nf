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

include { DROP  } from './workflows/drop'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_drop_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_drop_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_drop_pipeline'

include { samplesheetToList       } from 'plugin/nf-schema'

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

    def fasta = params.fasta ? Channel.value([id: 'fasta'], file(params.fasta)) : [[:], []]
    def fai   = params.fai   ? Channel.value([id: 'fai'], file(params.fai))     : [[:], []]
    // TODO accomodate for the ncbi and ucsc fasta

    def hpo_file = params.hpo_file ? Channel.value([id: 'hpo'], file(params.hpo_file)) : [[:], []]

    def gene_annotation = params.gene_annotation
        ? samplesheetToList(params.gene_annotation).collectEntries { name, gtf -> [ name, file(gtf) ] }
        : [:]

    def ec_gene_annotations = params.ec_gene_annotations ? params.ec_gene_annotations.tokenize(",") : []
    def ec_exclude_groups = params.ec_exclude_groups ? params.ec_exclude_groups.tokenize(",") : []

    def ae_groups = params.ae_groups ? params.ae_groups.tokenize(",") : []
    def ae_genes_to_test = params.ae_genes_to_test ? Utils.readYamlFile(file(params.ae_genes_to_test)) : [:]
    //
    // WORKFLOW: Run main workflow
    //
    DROP (
        // Global parameters
        PIPELINE_INITIALISATION.out.samplesheet, // derived from --input
        params.project_title,
        fasta,
        fai,
        gene_annotation,
        hpo_file,

        // Export counts parameters
        ec_gene_annotations,
        ec_exclude_groups,

        // Aberrant expression parameters
        params.ae_run,
        ae_groups,
        ae_genes_to_test,
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
