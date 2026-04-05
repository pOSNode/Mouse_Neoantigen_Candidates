process NETMHCPAN {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/netmhcpan/${meta.id}", mode: 'copy'

    // NetMHCpan-4.1 requires an academic licence from DTU Health Tech.
    // Build the image locally with docker/Dockerfile.netmhcpan after
    // downloading the binary from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
    container 'owais/netmhcpan:4.1'

    input:
    tuple val(meta), path(annotated_vcf), path(annotated_tbi)
    val   mhc_alleles   // comma-separated, e.g. "H-2-Kb,H-2-Db"

    output:
    tuple val(meta), path("*_netmhcpan_raw.tsv"),        emit: raw
    tuple val(meta), path("*_netmhcpan_candidates.tsv"), emit: candidates
    tuple val(meta), path("*_netmhcpan_summary.tsv"),    emit: summary
    path "versions.yml",                                  emit: versions

    script:
    def prefix   = meta.id
    def alleles  = mhc_alleles
    def lengths  = params.peptide_lengths
    def el_rank  = params.el_rank_cutoff
    def ba_rank  = params.ba_rank_cutoff
    """
    # 1. Extract mutant peptides from VEP VCF → FASTA
    python3 ${projectDir}/scripts/vcf_to_peptides.py \\
        --vcf   ${annotated_vcf} \\
        --out   ${prefix}_peptides.fa \\
        --flank ${params.peptide_flank}

    # 2. Run NetMHCpan-4.1 for each allele × length combination
    python3 ${projectDir}/scripts/run_netmhcpan.py \\
        --fasta    ${prefix}_peptides.fa \\
        --alleles  "${alleles}" \\
        --lengths  "${lengths}" \\
        --el-rank  ${el_rank} \\
        --ba-rank  ${ba_rank} \\
        --sample   ${prefix} \\
        --outdir   .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhcpan: \$(netMHCpan -version 2>&1 | grep "NetMHCpan" | head -1 | awk '{print \$2}')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
