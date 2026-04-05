#!/usr/bin/env python3
"""
neoantigen_predict.py
---------------------
Neoantigen Prediction Script

Steps:
  1. Parse VEP-annotated VCF → extract missense/frameshift variants
  2. Fetch mutant peptide sequences via Ensembl REST API
  3. Submit peptides to IEDB MHC-I and MHC-II binding prediction
  4. Filter by EL Rank threshold and write results

Author: Owais Siddiqi
"""

import argparse
import csv
import json
import logging
import sys
import time
from pathlib import Path

import pandas as pd
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"
IEDB_MHC_I   = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
IEDB_MHC_II  = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"

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


# ─── VCF Parsing ─────────────────────────────────────────────────────────────

def parse_vep_vcf(vcf_path: str) -> pd.DataFrame:
    """Extract VEP CSQ fields from an annotated VCF into a DataFrame."""
    records = []
    csq_fields = []

    with (
        open(vcf_path)
        if not vcf_path.endswith(".gz")
        else __import__("gzip").open(vcf_path, "rt")
    ) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                # Parse CSQ field names from header
                desc = line.split("Format: ")[1].rstrip('"\n>')
                csq_fields = desc.split("|")
                continue
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            chrom, pos, _, ref, alt, qual, filt, info = cols[:8]

            if filt not in ("PASS", "."):
                continue

            info_dict = {}
            for entry in info.split(";"):
                if "=" in entry:
                    k, v = entry.split("=", 1)
                    info_dict[k] = v

            csq_raw = info_dict.get("CSQ", "")
            for transcript in csq_raw.split(","):
                fields = transcript.split("|")
                if len(fields) != len(csq_fields):
                    continue
                csq = dict(zip(csq_fields, fields))

                consequences = set(csq.get("Consequence", "").split("&"))
                if not consequences & CONSEQUENCE_WHITELIST:
                    continue

                records.append(
                    {
                        "chrom":        chrom,
                        "pos":          int(pos),
                        "ref":          ref,
                        "alt":          alt,
                        "qual":         qual,
                        "gene":         csq.get("SYMBOL", ""),
                        "transcript":   csq.get("Feature", ""),
                        "protein_id":   csq.get("ENSP", ""),
                        "consequence":  csq.get("Consequence", ""),
                        "hgvsp":        csq.get("HGVSp", ""),
                        "aa_pos":       csq.get("Protein_position", ""),
                        "aa_change":    csq.get("Amino_acids", ""),
                    }
                )

    df = pd.DataFrame(records).drop_duplicates(subset=["protein_id", "hgvsp"])
    log.info(f"Parsed {len(df)} qualifying variants from VCF")
    return df


# ─── Peptide Fetching ─────────────────────────────────────────────────────────

def fetch_mutant_peptide(
    protein_id: str, aa_pos: str, aa_change: str, flank: int = 13
) -> str | None:
    """
    Fetch the canonical protein sequence from Ensembl REST API,
    apply the amino acid substitution and return a flanked peptide window.
    """
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
    try:
        pos = int(aa_pos.split("-")[0]) - 1  # convert to 0-indexed
    except ValueError:
        return None

    if pos >= len(seq) or seq[pos].upper() != ref_aa.upper():
        log.debug(f"AA mismatch at {protein_id}:{pos} — expected {ref_aa}, got {seq[pos]}")

    mut_seq = seq[:pos] + alt_aa + seq[pos + 1 :]
    start   = max(0, pos - flank)
    end     = min(len(mut_seq), pos + flank + 1)
    peptide = mut_seq[start:end]
    time.sleep(0.1)  # Ensembl rate-limit courtesy
    return peptide


# ─── IEDB Submission ──────────────────────────────────────────────────────────

def predict_mhc_i(peptide: str, alleles: list[str], lengths: list[int]) -> list[dict]:
    """Submit a peptide to IEDB MHC-I prediction endpoint."""
    results = []
    try:
        payload = {
            "method":        "recommended",
            "sequence_text": peptide,
            "allele":        ",".join(alleles),
            "length":        ",".join(map(str, lengths)),
        }
        r = requests.post(IEDB_MHC_I, data=payload, timeout=60)
        r.raise_for_status()
        reader = csv.DictReader(r.text.splitlines(), delimiter="\t")
        for row in reader:
            results.append(
                {
                    "allele":   row.get("allele", ""),
                    "peptide":  row.get("peptide", ""),
                    "core":     row.get("core", ""),
                    "el_score": float(row.get("ann_ic50", 0)),
                    "el_rank":  float(row.get("ann_rank", 999)),
                    "mhc_class": "I",
                }
            )
    except Exception as exc:
        log.warning(f"IEDB MHC-I prediction failed: {exc}")
    return results


def predict_mhc_ii(peptide: str, alleles: list[str]) -> list[dict]:
    """Submit a peptide to IEDB MHC-II prediction endpoint."""
    results = []
    try:
        payload = {
            "method":        "recommended",
            "sequence_text": peptide,
            "allele":        ",".join(alleles),
        }
        r = requests.post(IEDB_MHC_II, data=payload, timeout=60)
        r.raise_for_status()
        reader = csv.DictReader(r.text.splitlines(), delimiter="\t")
        for row in reader:
            results.append(
                {
                    "allele":   row.get("allele", ""),
                    "peptide":  row.get("sequence", ""),
                    "core":     row.get("core", ""),
                    "el_score": float(row.get("score", 0)),
                    "el_rank":  float(row.get("rank", 999)),
                    "mhc_class": "II",
                }
            )
    except Exception as exc:
        log.warning(f"IEDB MHC-II prediction failed: {exc}")
    return results


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Neoantigen prediction from VEP VCF")
    parser.add_argument("--vcf",            required=True,  help="VEP-annotated VCF (.vcf or .vcf.gz)")
    parser.add_argument("--sample",         required=True,  help="Sample ID")
    parser.add_argument("--alleles",        required=True,  help="Comma-separated MHC alleles")
    parser.add_argument("--lengths",        default="8,9,10,11", help="Peptide lengths for MHC-I")
    parser.add_argument("--el-rank-cutoff", type=float, default=2.0, help="%EL Rank cutoff")
    parser.add_argument("--outdir",         default=".",    help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_alleles = [a.strip() for a in args.alleles.split(",")]
    # Separate MHC-I vs MHC-II alleles (murine convention)
    mhc_i_alleles  = [a for a in all_alleles if "IA" not in a and "IE" not in a]
    mhc_ii_alleles = [a for a in all_alleles if "IA" in a or "IE" in a]
    lengths        = [int(l) for l in args.lengths.split(",")]

    # Parse VCF
    variants = parse_vep_vcf(args.vcf)
    if variants.empty:
        log.warning("No qualifying variants found. Writing empty outputs.")
        for suffix in ["_mhc_class_I_neoantigens.tsv",
                        "_mhc_class_II_neoantigens.tsv",
                        "_neoantigen_summary.tsv"]:
            (outdir / f"{args.sample}{suffix}").write_text("")
        return

    all_predictions: list[dict] = []

    for _, row in variants.iterrows():
        log.info(f"Processing {row['gene']} — {row['hgvsp']}")
        peptide = fetch_mutant_peptide(
            row["protein_id"], row["aa_pos"], row["aa_change"]
        )
        if not peptide:
            continue

        # MHC-I
        if mhc_i_alleles:
            preds = predict_mhc_i(peptide, mhc_i_alleles, lengths)
            for p in preds:
                p.update({
                    "gene":       row["gene"],
                    "hgvsp":      row["hgvsp"],
                    "consequence": row["consequence"],
                    "full_peptide_length": len(peptide),
                    "full_peptide": peptide,
                })
            all_predictions.extend(preds)

        # MHC-II
        if mhc_ii_alleles:
            preds = predict_mhc_ii(peptide, mhc_ii_alleles)
            for p in preds:
                p.update({
                    "gene":       row["gene"],
                    "hgvsp":      row["hgvsp"],
                    "consequence": row["consequence"],
                    "full_peptide_length": len(peptide),
                    "full_peptide": peptide,
                })
            all_predictions.extend(preds)

    df_all = pd.DataFrame(all_predictions)
    if df_all.empty:
        log.warning("No neoantigen predictions generated.")
        return

    # Filter by EL Rank
    df_pass = df_all[df_all["el_rank"] <= args.el_rank_cutoff].copy()
    df_pass.sort_values("el_rank", inplace=True)

    df_mhc_i  = df_pass[df_pass["mhc_class"] == "I"]
    df_mhc_ii = df_pass[df_pass["mhc_class"] == "II"]

    df_mhc_i.to_csv(outdir / f"{args.sample}_mhc_class_I_neoantigens.tsv",  sep="\t", index=False)
    df_mhc_ii.to_csv(outdir / f"{args.sample}_mhc_class_II_neoantigens.tsv", sep="\t", index=False)

    # Summary
    summary = pd.DataFrame(
        {
            "sample":          [args.sample],
            "total_variants":  [len(variants)],
            "peptides_tested": [df_all["hgvsp"].nunique()],
            "mhc_i_pass":      [len(df_mhc_i)],
            "mhc_ii_pass":     [len(df_mhc_ii)],
            "el_rank_cutoff":  [args.el_rank_cutoff],
            "top_mhc_i_gene":  [df_mhc_i["gene"].iloc[0] if not df_mhc_i.empty else "N/A"],
            "top_mhc_ii_gene": [df_mhc_ii["gene"].iloc[0] if not df_mhc_ii.empty else "N/A"],
        }
    )
    summary.to_csv(outdir / f"{args.sample}_neoantigen_summary.tsv", sep="\t", index=False)

    log.info(
        f"Done. MHC-I: {len(df_mhc_i)} candidates | MHC-II: {len(df_mhc_ii)} candidates"
    )


if __name__ == "__main__":
    main()
