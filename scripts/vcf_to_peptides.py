#!/usr/bin/env python3
"""
vcf_to_peptides.py
------------------
Parse a VEP-annotated VCF, fetch mutant protein sequences from the
Ensembl REST API, and write a FASTA file for NetMHCpan input.

Each FASTA record header encodes metadata needed downstream:
  >{sample}|{gene}|{transcript}|{hgvsp}|{consequence}

Author: Owais Siddiqi
"""

import argparse
import gzip
import logging
import sys
import time
from pathlib import Path

import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"

CONSEQUENCE_WHITELIST = {
    "missense_variant",
    "frameshift_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
}


def open_vcf(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_vep_vcf(vcf_path: str) -> list[dict]:
    records = []
    csq_fields = []

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                desc = line.split("Format: ")[1].rstrip('"\n>')
                csq_fields = desc.split("|")
                continue
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            chrom, pos, _, ref, alt, qual, filt, info = cols[:8]

            if filt not in ("PASS", "."):
                continue

            info_dict = {
                k: v
                for entry in info.split(";")
                if "=" in entry
                for k, v in [entry.split("=", 1)]
            }

            for transcript in info_dict.get("CSQ", "").split(","):
                fields = transcript.split("|")
                if len(fields) != len(csq_fields):
                    continue
                csq = dict(zip(csq_fields, fields))

                consequences = set(csq.get("Consequence", "").split("&"))
                if not consequences & CONSEQUENCE_WHITELIST:
                    continue

                records.append(
                    {
                        "chrom":       chrom,
                        "pos":         pos,
                        "ref":         ref,
                        "alt":         alt,
                        "gene":        csq.get("SYMBOL", ""),
                        "transcript":  csq.get("Feature", ""),
                        "protein_id":  csq.get("ENSP", ""),
                        "consequence": csq.get("Consequence", ""),
                        "hgvsp":       csq.get("HGVSp", ""),
                        "aa_pos":      csq.get("Protein_position", ""),
                        "aa_change":   csq.get("Amino_acids", ""),
                    }
                )

    # Deduplicate on protein_id + hgvsp
    seen = set()
    unique = []
    for r in records:
        key = (r["protein_id"], r["hgvsp"])
        if key not in seen and r["protein_id"]:
            seen.add(key)
            unique.append(r)

    log.info(f"Parsed {len(unique)} unique qualifying variants")
    return unique


def fetch_mutant_sequence(
    protein_id: str, aa_pos: str, aa_change: str, flank: int
) -> str | None:
    """Fetch canonical sequence from Ensembl and apply the mutation."""
    try:
        r = requests.get(
            f"{ENSEMBL_REST}/sequence/id/{protein_id}",
            params={"type": "protein", "content-type": "application/json"},
            timeout=30,
        )
        r.raise_for_status()
        seq = r.json().get("seq", "")
    except Exception as exc:
        log.warning(f"Ensembl fetch failed for {protein_id}: {exc}")
        return None

    if not seq or "/" not in aa_change:
        return None

    ref_aa, alt_aa = aa_change.split("/", 1)
    # Handle ranges e.g. "123-125"
    try:
        pos = int(aa_pos.split("-")[0]) - 1
    except ValueError:
        return None

    if pos >= len(seq):
        return None

    mut_seq = seq[:pos] + alt_aa + seq[pos + 1:]
    start   = max(0, pos - flank)
    end     = min(len(mut_seq), pos + flank + 1)
    peptide = mut_seq[start:end]

    time.sleep(0.1)  # Ensembl rate-limit courtesy
    return peptide


def main():
    parser = argparse.ArgumentParser(
        description="Extract mutant peptides from VEP VCF → FASTA for NetMHCpan"
    )
    parser.add_argument("--vcf",   required=True, help="VEP-annotated VCF (.vcf or .vcf.gz)")
    parser.add_argument("--out",   required=True, help="Output FASTA path")
    parser.add_argument("--flank", type=int, default=13,
                        help="AA flank around mutation (default: 13 → 27-mer)")
    args = parser.parse_args()

    variants = parse_vep_vcf(args.vcf)
    if not variants:
        log.warning("No qualifying variants — writing empty FASTA")
        Path(args.out).write_text("")
        return

    written = 0
    with open(args.out, "w") as fa:
        for v in variants:
            peptide = fetch_mutant_sequence(
                v["protein_id"], v["aa_pos"], v["aa_change"], args.flank
            )
            if not peptide:
                continue

            # Header encodes all metadata for downstream joining
            header = (
                f">{v['gene']}|{v['transcript']}|"
                f"{v['hgvsp']}|{v['consequence']}"
            )
            fa.write(f"{header}\n{peptide}\n")
            written += 1

    log.info(f"Written {written} peptide sequences to {args.out}")


if __name__ == "__main__":
    main()
