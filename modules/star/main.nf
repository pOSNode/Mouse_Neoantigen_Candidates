process STAR_INDEX {
    tag "GRCm39"
    label 'process_high_memory'
    publishDir "${params.outdir}/star_index", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.11a--h0033a41_0'

    input:
    path genome
    path gtf

    output:
    path "star_index/", emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index/ \\
        --genomeFastaFiles ${genome} \\
        --sjdbGTFfile ${gtf} \\
        --genomeSAindexNbases 14

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}

process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/star/${meta.id}", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.11a--h0033a41_0'

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam"),     emit: bam
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam.bai"), emit: bai
    tuple val(meta), path("*Log.final.out"),                      emit: log_final
    tuple val(meta), path("*SJ.out.tab"),                         emit: sj
    path "versions.yml",                                          emit: versions

    script:
    def prefix = meta.id
    """
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD \\
        --outSAMstrandField intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonical \\
        --outFileNamePrefix ${prefix}_ \\
        --outBAMsortingThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --limitBAMsortRAM ${task.memory.toBytes()} \\
        --outSAMattrRGline ID:${prefix} SM:${prefix} PL:ILLUMINA LB:${prefix}

    samtools index ${prefix}_Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
