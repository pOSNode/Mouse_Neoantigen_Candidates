process TRIMMOMATIC {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/trimmomatic/${meta.id}", mode: 'copy'

    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_paired_R{1,2}.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("*_unpaired_R{1,2}.fastq.gz"), emit: unpaired_reads
    tuple val(meta), path("*.log"),                       emit: log
    path "versions.yml",                                  emit: versions

    script:
    def prefix = meta.id
    """
    trimmomatic PE \\
        -threads ${task.cpus} \\
        -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${prefix}_paired_R1.fastq.gz   ${prefix}_unpaired_R1.fastq.gz \\
        ${prefix}_paired_R2.fastq.gz   ${prefix}_unpaired_R2.fastq.gz \\
        ILLUMINACLIP:${params.trim_adapters}:2:30:10:8:keepBothReads \\
        LEADING:${params.trim_leading} \\
        TRAILING:${params.trim_trailing} \\
        SLIDINGWINDOW:${params.trim_slidingwindow} \\
        MINLEN:${params.trim_minlen} \\
        2> ${prefix}_trimmomatic.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
