process PARABRICKS_FQ2BAM {
    tag "${meta.id}"
    label 'process_gpu'
    publishDir "${params.outdir}/parabricks/${meta.id}", mode: 'copy'

    // Requires NVIDIA Clara Parabricks licence — container pulled from NGC
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1'

    input:
    tuple val(meta), path(reads)
    path  star_index
    path  genome

    output:
    tuple val(meta), path("*.bam"),      emit: bam
    tuple val(meta), path("*.bam.bai"),  emit: bai
    tuple val(meta), path("*_duplicate_metrics.txt"), emit: metrics
    path "versions.yml",                 emit: versions

    script:
    def prefix = meta.id
    """
    pbrun rna_fq2bam \\
        --ref ${genome} \\
        --genome-lib-dir ${star_index} \\
        --in-fq ${reads[0]} ${reads[1]} \\
        --out-bam ${prefix}.bam \\
        --read-files-command zcat \\
        --two-pass-mode Basic \\
        --out-duplicate-metrics ${prefix}_duplicate_metrics.txt \\
        --num-threads ${task.cpus} \\
        --tmp-dir ./tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun version 2>&1 | grep "^Parabricks" | sed 's/Parabricks v//')
    END_VERSIONS
    """
}
