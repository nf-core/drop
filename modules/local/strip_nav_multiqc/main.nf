process STRIP_NAV_MULTIQC {
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e48219bad31e50a9c692218ddce90b8df9314f5844db308b05b6ada1fe9c6a54/data'
        : 'community.wave.seqera.io/library/beautifulsoup4_lxml_python:50fd3d504309f7fd'}"

    input:
    path(html)

    output:
    path ('*_mqc.html'), emit: report
    path ('versions.yml'), emit: versions

    script:
    def out = "${html.baseName}_mqc.html"

    """
    strip_nav.py --html ${html} --output ${out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch test_mqc.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
