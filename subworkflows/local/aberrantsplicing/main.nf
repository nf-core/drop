workflow ABERRANTSPLICING {
    take:
    inputs         // queue channel: [ val(meta), path(bam), path(bai), path(vcf), path(tbi), path(gene_count), path(splice_count) ]
    fraser_version // string:        Fraser version to use for aberrant splicing analysis
    include_groups // list:          A list of groups to include in the aberrant splicing analysis
    genes_to_test  // value channel: [ val(meta), path(genes_to_test) ]

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
            def header = ["sampleID\tbamFile\tPAIRED_END\tSTRAND"]
            def lines = []
            def index = 0
            metas.each {
                lines += [
                    metas[index].id,
                    bams[index].name,
                    metas[index].paired_end ? "TRUE" : "FALSE",
                    metas[index].strand ?: ""
                ].join("\t")
                index += 1
            }
            return [group, header + lines.sort()]
        }
        .collectFile(newLine:true, storeDir:"${params.outdir}/processed_data/aberrant_splicing/annotations/") { group, lines ->
            [ "${group}.tsv", lines.join("\n") ]
        }


    emit:
    versions = ch_versions
}
