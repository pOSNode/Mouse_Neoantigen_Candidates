process STRELKA2_SOMATIC {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/strelka2/${meta.id}", mode: 'copy'

    container 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'

    input:
    tuple val(meta), path(bam), path(bai)
    path  genome
    path  genome_fai

    output:
    tuple val(meta), path("*.somatic_snvs.vcf.gz"),          emit: vcf
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi"),      emit: tbi
    tuple val(meta), path("*.somatic_indels.vcf.gz"),         emit: indels_vcf
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"),     emit: indels_tbi
    path "versions.yml",                                       emit: versions

    script:
    def prefix = meta.id
    """
    # Configure Strelka2 run
    configureStrelkaGermlineWorkflow.py \\
        --bam ${bam} \\
        --referenceFasta ${genome} \\
        --rna \\
        --runDir strelka_run

    # Execute workflow
    strelka_run/runWorkflow.py \\
        -m local \\
        -j ${task.cpus}

    # Rename outputs with sample prefix
    mv strelka_run/results/variants/genome.S1.vcf.gz      ${prefix}.somatic_snvs.vcf.gz
    mv strelka_run/results/variants/genome.S1.vcf.gz.tbi  ${prefix}.somatic_snvs.vcf.gz.tbi
    mv strelka_run/results/variants/genome.S1.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz  2>/dev/null || true
    mv strelka_run/results/variants/genome.S1.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$(configureStrelkaGermlineWorkflow.py --version)
    END_VERSIONS
    """
}
