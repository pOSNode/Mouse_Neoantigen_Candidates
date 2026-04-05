#!/usr/bin/env python3
"""
run_netmhcpan.py
----------------
NetMHCpan-4.1 Runner & Candidate Scorer

Workflow:
  1. For each allele × peptide-length combination, invoke netMHCpan
  2. Parse raw output → unified DataFrame
  3. Filter on %EL Rank (presentation) and %BA Rank (binding affinity)
  4. Score and rank candidates using a composite priority score
  5. Write:
       {sample}_netmhcpan_raw.tsv        — all predictions pre-filter
       {sample}_netmhcpan_candidates.tsv — PASS candidates, ranked
       {sample}_netmhcpan_summary.tsv    — per-gene/allele summary

Composite score (higher = better candidate):
    priority_score = (1 - EL_rank/100) * 0.6
                   + (1 - BA_rank/100) * 0.3
                   + (1 / peptide_length) * 0.1

NetMHCpan binary must be on PATH (netMHCpan).
Academic licence: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/

Author: Owais Siddiqi
"""

import argparse
import logging
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ─── NetMHCpan output column spec ────────────────────────────────────────────
# Columns emitted by netMHCpan -xls output (positional, tab-separated)
NMH_COLS = [
    "pos", "mhc", "peptide", "core", "of", "gp", "gl",
    "ip", "il", "icore", "identity", "score_el", "el_rank",
    "score_ba", "ba_rank", "ave", "nmer", "bind_level",
]


# ─── Runner ──────────────────────────────────────────────────────────────────

def run_netmhcpan(
    fasta: str,
    allele: str,
    length: int,
    tmp_dir: Path,
) -> pd.DataFrame:
    """
    Invoke netMHCpan for a single allele/length and parse the XLS output.
    Returns a DataFrame of raw predictions.
    """
    xls_out = tmp_dir / f"{allele}_{length}.xls"

    cmd = [
        "netMHCpan",
        "-f", fasta,
        "-a", allele,
        "-l", str(length),
        "-BA",          # include binding affinity predictions
        "-xls",
        "-xlsfile", str(xls_out),
    ]

    log.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        log.warning(
            f"netMHCpan failed for {allele} / {length}-mer:\n{exc.stderr}"
        )
        return pd.DataFrame()
    except FileNotFoundError:
        log.error("netMHCpan binary not found on PATH")
        raise

    if not xls_out.exists():
        log.warning(f"No XLS output produced for {allele} / {length}-mer")
        return pd.DataFrame()

    df = _parse_xls(xls_out, allele, length)
    log.info(f"  {allele} / {length}-mer → {len(df)} predictions")
    return df


def _parse_xls(xls_path: Path, allele: str, length: int) -> pd.DataFrame:
    """
    Parse netMHCpan tab-delimited XLS output into a clean DataFrame.
    Skips header/footer comment lines (starting with '#' or 'Pos').
    """
    rows = []
    with open(xls_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#") or line.startswith("Pos"):
                continue
            parts = line.split("\t")
            # Pad or truncate to expected column count
            parts += [""] * (len(NMH_COLS) - len(parts))
            rows.append(parts[: len(NMH_COLS)])

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows, columns=NMH_COLS)

    # Type coercions
    for col in ["el_rank", "ba_rank", "score_el", "score_ba"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["peptide_length"] = length
    df["allele"]         = allele  # overwrite abbreviated form with canonical

    # Parse identity → gene, transcript, hgvsp, consequence
    meta = df["identity"].str.split("|", expand=True)
    df["gene"]        = meta[0] if 0 in meta else ""
    df["transcript"]  = meta[1] if 1 in meta else ""
    df["hgvsp"]       = meta[2] if 2 in meta else ""
    df["consequence"] = meta[3] if 3 in meta else ""

    return df


# ─── Scoring ─────────────────────────────────────────────────────────────────

def compute_priority_score(df: pd.DataFrame) -> pd.Series:
    """
    Composite priority score combining EL rank, BA rank, and peptide length.
    Range: 0 → 1 (higher = better candidate for follow-up).
    """
    el  = 1 - (df["el_rank"].clip(0, 100) / 100)
    ba  = 1 - (df["ba_rank"].clip(0, 100) / 100)
    ln  = 1 / df["peptide_length"].clip(lower=1)
    return (el * 0.6 + ba * 0.3 + ln * 0.1).round(6)


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run NetMHCpan-4.1 and score neoantigen candidates"
    )
    parser.add_argument("--fasta",    required=True,  help="Mutant peptide FASTA")
    parser.add_argument("--alleles",  required=True,  help="Comma-separated MHC-I alleles")
    parser.add_argument("--lengths",  default="8,9,10,11", help="Peptide lengths")
    parser.add_argument("--el-rank",  type=float, default=2.0,
                        help="%EL Rank cutoff (default: 2.0)")
    parser.add_argument("--ba-rank",  type=float, default=2.0,
                        help="%BA Rank cutoff (default: 2.0, ~500 nM)")
    parser.add_argument("--sample",   required=True, help="Sample ID")
    parser.add_argument("--outdir",   default=".",   help="Output directory")
    args = parser.parse_args()

    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    alleles = [a.strip() for a in args.alleles.split(",")]
    lengths = [int(l) for l in args.lengths.split(",")]

    fasta = args.fasta
    if not Path(fasta).exists() or Path(fasta).stat().st_size == 0:
        log.warning("Empty FASTA — no peptides to score. Writing empty outputs.")
        for suffix in ["_netmhcpan_raw.tsv",
                        "_netmhcpan_candidates.tsv",
                        "_netmhcpan_summary.tsv"]:
            (outdir / f"{args.sample}{suffix}").write_text("")
        return

    # ── Run NetMHCpan for all allele × length combinations ──────────────────
    all_frames = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for allele in alleles:
            for length in lengths:
                df = run_netmhcpan(fasta, allele, length, tmp_dir)
                if not df.empty:
                    all_frames.append(df)

    if not all_frames:
        log.warning("No predictions returned from NetMHCpan.")
        return

    df_raw = pd.concat(all_frames, ignore_index=True)
    df_raw["sample"] = args.sample

    # ── Write raw ────────────────────────────────────────────────────────────
    raw_cols = [
        "sample", "gene", "hgvsp", "consequence",
        "allele", "peptide_length", "peptide", "core",
        "score_el", "el_rank", "score_ba", "ba_rank", "bind_level",
    ]
    df_raw[raw_cols].to_csv(
        outdir / f"{args.sample}_netmhcpan_raw.tsv", sep="\t", index=False
    )
    log.info(f"Raw predictions: {len(df_raw)}")

    # ── Filter candidates ────────────────────────────────────────────────────
    df_pass = df_raw[
        (df_raw["el_rank"] <= args.el_rank) &
        (df_raw["ba_rank"] <= args.ba_rank)
    ].copy()

    if df_pass.empty:
        log.warning(
            f"No candidates passed EL rank ≤ {args.el_rank} "
            f"AND BA rank ≤ {args.ba_rank}"
        )
    else:
        df_pass["priority_score"] = compute_priority_score(df_pass)
        df_pass.sort_values("priority_score", ascending=False, inplace=True)
        df_pass.insert(0, "rank", range(1, len(df_pass) + 1))

        candidate_cols = [
            "rank", "sample", "gene", "hgvsp", "consequence",
            "allele", "peptide_length", "peptide", "core",
            "el_rank", "ba_rank", "priority_score", "bind_level",
        ]
        df_pass[candidate_cols].to_csv(
            outdir / f"{args.sample}_netmhcpan_candidates.tsv",
            sep="\t", index=False,
        )
        log.info(f"Candidates (EL ≤ {args.el_rank}, BA ≤ {args.ba_rank}): {len(df_pass)}")

    # ── Summary ──────────────────────────────────────────────────────────────
    summary_rows = []
    for gene, gdf in df_pass.groupby("gene") if not df_pass.empty else []:
        for allele, adf in gdf.groupby("allele"):
            best = adf.sort_values("el_rank").iloc[0]
            summary_rows.append(
                {
                    "sample":          args.sample,
                    "gene":            gene,
                    "allele":          allele,
                    "n_candidates":    len(adf),
                    "best_peptide":    best["peptide"],
                    "best_el_rank":    best["el_rank"],
                    "best_ba_rank":    best["ba_rank"],
                    "best_priority":   best["priority_score"],
                    "hgvsp":           best["hgvsp"],
                }
            )

    pd.DataFrame(summary_rows).to_csv(
        outdir / f"{args.sample}_netmhcpan_summary.tsv", sep="\t", index=False
    )

    log.info("Done.")


if __name__ == "__main__":
    main()
