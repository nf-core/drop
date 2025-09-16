include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats'

workflow BAM_STATS_IDXSTATS_MERGE {
    take:
    bams            // queue channel:   [ val(meta), path(bam), path(bai) ]
    externals       // queue channel:    val(meta)
    include_groups  // list:            A list of groups to include in the bam stats


    main:
    def ucsc2ncbi = file("${projectDir}/assets/chr_UCSC_NCBI.txt", checkIfExists:true).splitCsv(header: false, sep: " ")
    def ucsc_chr = ucsc2ncbi.collect { it[0] }
    def ncbi_chr = ucsc2ncbi.collect { it[1] }
    def ch_versions = Channel.empty()

    SAMTOOLS_IDXSTATS(
        bams
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    def externals_entries = externals
        .map { meta -> [ meta, "${meta.id}\tNA" ] }

    SAMTOOLS_IDXSTATS.out.idxstats
        // Get the total read count for each sample
        .map { meta, idxstats ->
            def split_idxstats = idxstats.splitCsv(header: false, sep: "\t")
            // Check which format the chromosome names are in (UCSC = chrX, NCBI = X)
            def chr_names = split_idxstats[0][0].startsWith("chr") ? ucsc_chr : ncbi_chr
            def total_mapped_reads = 0
            split_idxstats.each { line ->
                // Only use the reads that are in the main chromosomes
                if (!chr_names.contains(line[0])) {
                    return
                }
                total_mapped_reads += line[2].toInteger()
            }
            return [ meta, "${meta.id}\t${total_mapped_reads}" ]
        }
        .tap { ch_individuals }
        .collectFile(newLine:true, storeDir:"${params.outdir}/processed_data/aberrant_expression/bam_stats/") { meta, line ->
            [ "${meta.id}.txt", line]
        }

    ch_merged_bam_stats = ch_individuals
        .mix(externals_entries)
        .map { meta, line ->
            def groups = meta.drop_group.tokenize(",").intersect(include_groups)
            [ groups, meta, line ]
        }
        .transpose(by: 0)
        .map { group, meta, line ->
            [ groupKey(group, meta.drop_group_counts.get(group)), line ]
        }
        .groupTuple()
        // Add the file header and sort the lines on sample name
        .map { group, lines ->
            return [ group, (["sampleID\trecord_count"] + lines.sort()).join("\n") ]
        }
        // Create the merged bam stats file
        .collectFile(newLine:true, storeDir:"${params.outdir}/processed_data/aberrant_expression/bam_stats/") { group, line ->
            [ "${group}.txt", line ]
        }
        .map {
            path -> [ path.baseName.replaceFirst(/\.txt$/,'') , path ]
        }

    emit:
    versions =         ch_versions          // value channel: path(versions)
    merged_bam_stats = ch_merged_bam_stats  // value channel: path(merged_bam_stats)
}
