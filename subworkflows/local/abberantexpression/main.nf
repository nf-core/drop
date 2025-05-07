workflow ABBERANTEXPRESSION {
    take:
    bams // queue channel: [ val(meta), path(bam), path(bai) ]

    main:
    def ch_versions = Channel.empty()
    // TODO write the subworkflow

    emit:
    versions = ch_versions
}
