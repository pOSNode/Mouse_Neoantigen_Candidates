#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    RNA-Seq Neoantigen Discovery Pipeline
========================================================================================
    FastQ → QC → Alignment (STAR/Parabricks) → BAM Post-Processing
    → Variant Calling (Strelka2) → Filtering (GATK) → Annotation (VEP)
    → Neoantigen Scoring (NetMHCpan-4.1)
    GRCm39 (Mus musculus)
----------------------------------------------------------------------------------------
*/

include { FASTQC          } from './modules/fastqc/main'
include { TRIMMOMATIC     } from './modules/trimmomatic/main'
include { STAR_INDEX      } from './modules/star/main'
include { STAR_ALIGN      } from './modules/star/main'
include { PARABRICKS_FQ2BAM } from './modules/parabricks/main'
include { STRELKA2_SOMATIC  } from './modules/strelka2/main'
include { GATK_FILTRATION   } from './modules/gatk/main'
include { VEP_ANNOTATE      } from './modules/vep/main'
include { NETMHCPAN          } from './modules/netmhcpan/main'

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {

    // ── Input channels ──────────────────────────────────────────────────────────
    ch_reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, files -> [ [id: sample_id], files ] }

    ch_genome    = file(params.genome,    checkIfExists: true)
    ch_genome_fai = file(params.genome_fai, checkIfExists: true)
    ch_genome_dict = file(params.genome_dict, checkIfExists: true)
    ch_star_index = file(params.star_index, checkIfExists: true)

    // ── QC ───────────────────────────────────────────────────────────────────────
    FASTQC(ch_reads)
    TRIMMOMATIC(ch_reads)

    // ── Alignment ────────────────────────────────────────────────────────────────
    if (params.use_parabricks) {
        PARABRICKS_FQ2BAM(
            TRIMMOMATIC.out.trimmed_reads,
            ch_star_index,
            ch_genome
        )
        ch_bam = PARABRICKS_FQ2BAM.out.bam
    } else {
        STAR_ALIGN(
            TRIMMOMATIC.out.trimmed_reads,
            ch_star_index
        )
        ch_bam = STAR_ALIGN.out.bam
    }

    // ── Variant Calling ──────────────────────────────────────────────────────────
    STRELKA2_SOMATIC(
        ch_bam,
        ch_genome,
        ch_genome_fai
    )

    // ── Filtering ────────────────────────────────────────────────────────────────
    GATK_FILTRATION(
        STRELKA2_SOMATIC.out.vcf,
        ch_genome,
        ch_genome_fai,
        ch_genome_dict
    )

    // ── Annotation ───────────────────────────────────────────────────────────────
    VEP_ANNOTATE(
        GATK_FILTRATION.out.filtered_vcf,
        ch_genome,
        ch_genome_fai
    )

    // ── Neoantigen Prediction + Scoring (NetMHCpan-4.1) ─────────────────────────
    NETMHCPAN(
        VEP_ANNOTATE.out.annotated_vcf,
        params.mhc_alleles
    )
}
