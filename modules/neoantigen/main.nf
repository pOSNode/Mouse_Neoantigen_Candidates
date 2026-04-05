process NEOANTIGEN_PREDICT {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/neoantigens/${meta.id}", mode: 'copy'

    container 'python:3.11-slim'

    input:
    tuple val(meta), path(annotated_vcf), path(annotated_tbi)
    val   mhc_alleles

    output:
    tuple val(meta), path("*_mhc_class_I_neoantigens.tsv"),  emit: mhc_i
    tuple val(meta), path("*_mhc_class_II_neoantigens.tsv"), emit: mhc_ii
    tuple val(meta), path("*_neoantigen_summary.tsv"),        emit: summary
    path "versions.yml",                                       emit: versions

    script:
    def prefix  = meta.id
    def alleles = mhc_alleles
    def lengths = params.peptide_lengths
    def el_rank = params.el_rank_cutoff
    """
    pip install requests pandas --quiet

    python ${projectDir}/scripts/neoantigen_predict.py \\
        --vcf ${annotated_vcf} \\
        --sample ${prefix} \\
        --alleles "${alleles}" \\
        --lengths "${lengths}" \\
        --el-rank-cutoff ${el_rank} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
