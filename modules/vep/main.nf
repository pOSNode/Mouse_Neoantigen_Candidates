process VEP_ANNOTATE {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/vep/${meta.id}", mode: 'copy'

    container 'ensemblorg/ensembl-vep:release_111.0'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  genome
    path  genome_fai

    output:
    tuple val(meta), path("*.vep.vcf.gz"),     emit: annotated_vcf
    tuple val(meta), path("*.vep.vcf.gz.tbi"), emit: annotated_tbi
    tuple val(meta), path("*.vep_summary.html"), emit: summary
    path "versions.yml",                         emit: versions

    script:
    def prefix = meta.id
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${prefix}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --species mus_musculus \\
        --assembly GRCm39 \\
        --cache \\
        --offline \\
        --fasta ${genome} \\
        --fork ${task.cpus} \\
        --everything \\
        --canonical \\
        --hgvs \\
        --protein \\
        --uniprot \\
        --domains \\
        --af_gnomade \\
        --stats_file ${prefix}.vep_summary.html

    bgzip -c ${prefix}.vep.vcf > ${prefix}.vep.vcf.gz
    tabix -p vcf ${prefix}.vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$(echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
