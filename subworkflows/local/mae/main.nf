include { MAE_CREATESNVS } from '../../../modules/local/mae/createsnvs/main'

workflow MAE {
    take:
    input           // queue channel: [ val(meta), path(vcf), path(tbi), path(bam), path(bai) ]
    include_groups  // list         : A list of groups to include in the mono allelic expression analysis
    ncbi2ucsc       // value channel: path to the NCBI to UCSC mapping file
    ucsc2ncbi       // value channel: path to the UCSC to NCBI mapping file

    main:
    def ch_versions = Channel.empty()

    def ch_filtered_inputs = Channel.empty()
    if (include_groups) {
        ch_filtered_inputs = input.filter { meta, _vcf, _tbi, _bam, _bai ->
            meta.drop_group.tokenize(",").intersect(include_groups).size() > 0
        }
    } else {
        ch_filtered_inputs = input
    }

    def ch_createsnvs_input = ch_filtered_inputs
        .filter { _meta, vcf, tbi, bam, bai ->
            vcf && tbi && bam && bai
        }
        .map { meta, vcf, tbi, bam, bai ->
            [ meta, vcf, tbi, bam, bai, meta.id ]
        }

    MAE_CREATESNVS(
        ch_createsnvs_input,
        ncbi2ucsc,
        ucsc2ncbi
    )
    ch_versions = ch_versions.mix(MAE_CREATESNVS.out.versions.first())

    emit:
    versions = ch_versions
}
