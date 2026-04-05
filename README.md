# RNA-Seq Neoantigen Discovery Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![GRCm39](https://img.shields.io/badge/genome-GRCm39-4a90e2)](https://www.ensembl.org/Mus_musculus)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0006--4754--8200-a6ce39?logo=orcid)](https://orcid.org/0009-0006-4754-8200)

End-to-end Nextflow DSL2 pipeline for somatic variant calling and MHC-I neoantigen prediction from murine RNA-Seq data (GRCm39 / *Mus musculus*). Designed for personalised cancer immunotherapy research.

---

## Pipeline Overview

```
FastQ → FastQC + Trimmomatic → STAR (or NVIDIA Parabricks)
     → Strelka2 (Somatic/Germline) → GATK Hard-Filtering
     → VEP Annotation → Mutant Peptide Extraction
     → NetMHCpan-4.1 Scoring → Ranked Neoantigen Candidates
```

| Stage | Tool | Output |
|---|---|---|
| Quality Control | FastQC, Trimmomatic | QC reports, trimmed FASTQ |
| Alignment | STAR 2.7.11 / NVIDIA Parabricks `rna_fq2bam` | Sorted, indexed BAM |
| Variant Calling | Strelka2 2.9.10 | Somatic + germline VCF |
| Variant Filtering | GATK 4.5 `VariantFiltration` | PASS-only VCF |
| Annotation | Ensembl VEP 111 | Functionally annotated VCF |
| Peptide Extraction | Ensembl REST API + `vcf_to_peptides.py` | Mutant peptide FASTA |
| Neoantigen Scoring | NetMHCpan-4.1 + `run_netmhcpan.py` | Ranked MHC-I candidate TSVs |

---

## Requirements

- **Nextflow** >= 23.04.0
- **Docker** or **Singularity**
- **Java** >= 17
- **NetMHCpan-4.1** — academic licence required (see [Docker setup](#docker-setup))
- *(Optional)* NVIDIA GPU + Clara Parabricks licence for GPU-accelerated alignment

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/pOSNode/rna-neoantigen-pipeline.git
cd rna-neoantigen-pipeline

# 2. Run test profile (uses bundled chr19 test data)
nextflow run main.nf -profile test,docker

# 3. Full run
nextflow run main.nf \
    -profile docker \
    --reads  '/data/fastq/*_R{1,2}.fastq.gz' \
    --outdir results
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--reads` | `null` | Glob pattern for paired FASTQ files |
| `--outdir` | `results` | Output directory |
| `--genome` | — | Reference FASTA (GRCm39) |
| `--star_index` | — | Pre-built STAR index directory |
| `--use_parabricks` | `false` | Use GPU-accelerated alignment |
| `--qual_threshold` | `30` | GATK QUAL hard-filter threshold |
| `--dp_threshold` | `10` | GATK DP hard-filter threshold |
| `--mhc_alleles` | `H-2-Kb,H-2-Db` | Murine MHC-I alleles |
| `--peptide_lengths` | `8,9,10,11` | MHC-I peptide lengths |
| `--peptide_flank` | `13` | AA flank around mutation site |
| `--el_rank_cutoff` | `2.0` | %EL Rank filter (<=0.5 strong binder) |
| `--ba_rank_cutoff` | `2.0` | %BA Rank filter (~500 nM IC50) |

---

## Profiles

| Profile | Description |
|---|---|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers (HPC) |
| `gcp` | Google Cloud Batch executor |
| `hpc` | SLURM executor |
| `test` | Local test run with bundled chr19 data |

---

## Docker Setup

All containers except NetMHCpan pull automatically from public registries when running with `-profile docker` or `-profile singularity`.

### NetMHCpan-4.1 — manual build required

> **Licence notice:** NetMHCpan-4.1 is free for academic, non-commercial use only. The binary is **not** redistributable and therefore **not** included in this repository and **not** hosted on any public container registry. You must obtain your own copy directly from DTU Health Tech and build the image locally before running the neoantigen scoring step.

```bash
# 1. Register and download the Linux binary (free for academic use):
#    https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
#    Place the downloaded tarball at: docker/netMHCpan-4.1b.Linux.tar.gz

# 2. Build the image locally — this must be done before running the pipeline
docker build -f docker/Dockerfile.netmhcpan -t owais/netmhcpan:4.1 .

# 3. Verify
docker run --rm owais/netmhcpan:4.1 netMHCpan -version
```

If you attempt to run the pipeline without building this image first, the `NETMHCPAN` process will fail with a container-not-found error.

---

## Outputs

```
results/
├── fastqc/               # Per-sample FastQC HTML + ZIP reports
├── trimmomatic/          # Trimmed FASTQs + trimming logs
├── star/                 # Aligned BAMs + STAR logs
├── strelka2/             # Somatic + germline VCFs
├── gatk/                 # PASS-filtered VCFs
├── vep/                  # VEP-annotated VCFs + HTML summaries
├── netmhcpan/
│   ├── *_netmhcpan_raw.tsv             # All predictions (pre-filter)
│   ├── *_netmhcpan_candidates.tsv      # PASS candidates, ranked by priority score
│   └── *_netmhcpan_summary.tsv         # Per-gene/allele best-binder summary
└── pipeline_info/        # Nextflow timeline, report, trace, DAG
```

### Candidate Output Columns

| Column | Description |
|---|---|
| `rank` | Priority rank across all allele/length combinations |
| `gene` | Gene symbol (e.g. `Hsp90ab1`, `Ncl`) |
| `allele` | MHC-I allele (e.g. `H-2-Kb`) |
| `peptide_length` | Mer length (8-11) |
| `peptide` | Presented peptide sequence |
| `core` | NetMHCpan binding core (9-mer) |
| `el_rank` | %EL Rank — predicted ligand likelihood (<=0.5 strong binder) |
| `ba_rank` | %BA Rank — predicted binding affinity rank (<=2.0 ~500 nM) |
| `priority_score` | Composite score: 0.6x(1-EL/100) + 0.3x(1-BA/100) + 0.1x(1/len) |
| `hgvsp` | HGVSp notation of the causal variant |
| `consequence` | VEP consequence term |
| `bind_level` | NetMHCpan classification: SB strong / WB weak binder |

---

## Citation

If you use this pipeline, please cite:

> Siddiqi O. (2024). *RNA-Seq Neoantigen Discovery Pipeline* (v1.0.0). GitHub. https://github.com/pOSNode/rna-neoantigen-pipeline
> ORCID: 0009-0006-4754-8200

---

## Licence

MIT © Owais Siddiqi
