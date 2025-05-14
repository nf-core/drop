include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats'

workflow BAM_STATS_IDXSTATS {
    take:
    bams        // queue channel:   [ val(meta), path(bam), path(bai) ]
    exclude_ids // value channel:   A list of ids to exclude from the bam stats


    main:
    def ucsc2ncbi = file("${projectDir}/assets/chr_UCSC_NCBI.txt", checkIfExists:true).splitCsv(header: false, sep: " ")
    def ucsc_chr = ucsc2ncbi.collect { it[0] }
    def ncbi_chr = ucsc2ncbi.collect { it[1] }
    def ch_versions = Channel.empty()

    SAMTOOLS_IDXSTATS(
        bams
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    def exclude_ids_entries = exclude_ids
        .flatten()
        .map { id ->
            return "${id}\tNA"
        }

    def ch_merged_bam_stats = SAMTOOLS_IDXSTATS.out.idxstats
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
            return "${meta.id}\t${total_mapped_reads}"
        }
        // Mix in the excluded IDs that will have unknown read counts
        .mix(exclude_ids_entries)
        .collect()
        // Add the file header and sort the lines on sample name
        .map { lines ->
            return ["sampleID\trecord_count"] + lines.sort()
        }
        .flatten()
        // Create the merged bam stats file
        .collectFile(name: "merged_bam_stats.tsv", newLine:true, sort:false) // TODO this should be another name + add storeDir with the output dir location
        .collect()


    emit:
    versions =         ch_versions          // value channel: path(versions)
    merged_bam_stats = ch_merged_bam_stats  // value channel: path(merged_bam_stats)
}
