process GATK_FILTRATION {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/gatk/${meta.id}", mode: 'copy'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  genome
    path  genome_fai
    path  genome_dict

    output:
    tuple val(meta), path("*.hard_filtered.vcf.gz"),     emit: filtered_vcf
    tuple val(meta), path("*.hard_filtered.vcf.gz.tbi"), emit: filtered_tbi
    path "versions.yml",                                  emit: versions

    script:
    def prefix = meta.id
    def qual    = params.qual_threshold
    def dp      = params.dp_threshold
    """
    # Decompress for GATK input
    gatk VariantFiltration \\
        -V ${vcf} \\
        -R ${genome} \\
        --filter-expression "QUAL < ${qual}" --filter-name "LowQual" \\
        --filter-expression "DP < ${dp}"   --filter-name "LowDepth" \\
        -O ${prefix}.hard_filtered_raw.vcf.gz

    # Remove flagged variants (keep only PASS)
    gatk SelectVariants \\
        -V ${prefix}.hard_filtered_raw.vcf.gz \\
        --exclude-filtered \\
        -O ${prefix}.hard_filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep "^GATK" | sed 's/GATK v//')
    END_VERSIONS
    """
}
