process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc/${meta.id}", mode: 'copy'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path "versions.yml",             emit: versions

    script:
    def prefix = meta.id
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
