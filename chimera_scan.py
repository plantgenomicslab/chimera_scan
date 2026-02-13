#!/usr/bin/env python3
"""Standalone whole-proteome chimera scan pipeline.

Scans an entire reference proteome for chimeric (fused) genes by:
  1. MMseqs2 search: reference proteins vs SwissProt
  2. Streaming Pattern B scan (>80% coverage filter): bimodal target detection
  3. Self-blast: candidates vs reference proteome (self-excluded), Pattern B
  4. SwissProt description annotation with N/C-term start/end + distinct gene count

Usage:
  python chimera_scan.py \
      --ref-fasta reference_proteins.fasta \
      --swissprot-db ../../mmseqs_db/uniprot_sprot_odb10_plants \
      --swissprot-fasta swissprot_plants.fasta \
      --mmseqs-path mmseqs --mmseqs-threads 16 --verbose
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MMSEQS_COLUMNS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
    "pident", "qcov", "tcov", "qlen", "tlen",
]

MMSEQS_DTYPES = {
    "query": "str", "target": "str",
    "fident": "float32", "alnlen": "int32", "mismatch": "int32",
    "gapopen": "int32", "qstart": "int32", "qend": "int32",
    "tstart": "int32", "tend": "int32", "evalue": "float64",
    "bits": "float32", "pident": "float32", "qcov": "float32",
    "tcov": "float32", "qlen": "int32", "tlen": "int32",
}

COLS_HIT = ["query", "target", "evalue", "alnlen", "qstart", "qend",
            "tstart", "tend", "qlen", "tlen", "fident", "bits"]

# Pattern B thresholds
NTERM_CUTOFF = 0.35
CTERM_CUTOFF = 0.65
MIN_DOMAIN_TARGETS = 3

SWISSPROT_COLS = ["query", "target", "fident", "alnlen", "qstart", "qend",
                  "tstart", "tend", "evalue", "bits", "qlen", "tlen", "theader"]


# ---------------------------------------------------------------------------
# Interval utilities
# ---------------------------------------------------------------------------

def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    sorted_iv = sorted(intervals)
    merged = [list(sorted_iv[0])]
    for s, e in sorted_iv[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(m) for m in merged]


def intervals_coverage(intervals: List[Tuple[int, int]], seqlen: int) -> float:
    if seqlen <= 0:
        return 0.0
    return min(sum(e - s for s, e in intervals) / seqlen, 1.0)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def detect_header(path: str) -> bool:
    with open(path) as fh:
        return fh.readline().strip().startswith("query")


def extract_fasta_subset(in_fasta: str, ids: Set[str],
                         out_fasta: str) -> int:
    """Extract sequences matching ids from a FASTA file. Returns count."""
    from Bio import SeqIO
    count = 0
    with open(out_fasta, "w") as fh:
        for rec in SeqIO.parse(in_fasta, "fasta"):
            if rec.id in ids:
                SeqIO.write(rec, fh, "fasta")
                count += 1
    return count


def read_mmseqs_result(path: str, evalue_cutoff: float,
                       query_filter: Optional[Set[str]] = None,
                       exclude_self: bool = False) -> pd.DataFrame:
    """Read an MMseqs2 result file, applying filters."""
    has_hdr = detect_header(path)
    df = pd.read_csv(
        path, sep="\t",
        header=0 if has_hdr else None,
        names=None if has_hdr else MMSEQS_COLUMNS,
        usecols=COLS_HIT,
        dtype={c: MMSEQS_DTYPES[c] for c in COLS_HIT},
    )
    df = df[df["evalue"] <= evalue_cutoff]
    if query_filter is not None:
        df = df[df["query"].isin(query_filter)]
    if exclude_self:
        df = df[df["query"] != df["target"]]
    return df


# ---------------------------------------------------------------------------
# MMseqs2 runner
# ---------------------------------------------------------------------------

def run_mmseqs_search(query_fasta: str, target: str,
                      output_tsv: str, mmseqs_path: str,
                      threads: int, tmp_dir: str,
                      evalue: float = 1e-3,
                      verbose: bool = False,
                      fmt: Optional[str] = None) -> str:
    """Run mmseqs easy-search. Reuses existing result if present."""
    if os.path.isfile(output_tsv) and os.path.getsize(output_tsv) > 0:
        print(f"  Reusing existing result: {output_tsv}")
        return output_tsv

    if fmt is None:
        fmt = ",".join(MMSEQS_COLUMNS)
    cmd = [
        mmseqs_path, "easy-search",
        query_fasta, target, output_tsv, tmp_dir,
        "--format-output", fmt,
        "--threads", str(threads),
        "-e", str(evalue),
        "--max-seqs", "500",
    ]
    print(f"  Running: {' '.join(cmd[:6])} ...", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  MMseqs2 STDERR:\n{result.stderr}", file=sys.stderr)
        sys.exit(f"MMseqs2 failed with exit code {result.returncode}")
    if verbose and result.stderr:
        for line in result.stderr.strip().split("\n")[-3:]:
            print(f"    {line}")
    return output_tsv


# ---------------------------------------------------------------------------
# Pattern B: single query analysis
# ---------------------------------------------------------------------------

def _analyze_single_query_b(query_id: str, hits: pd.DataFrame,
                             min_hits: int
                             ) -> Optional[Tuple[dict, pd.DataFrame]]:
    """Test a single query for bimodal target distribution (chimera signal).

    Returns (summary_dict, annotated_hits_df) or None.
    """
    qlen = int(hits["qlen"].iloc[0])
    if len(hits) < min_hits:
        return None

    hits_c = hits.copy()
    hits_c["qmid"] = (hits_c["qstart"] + hits_c["qend"]) / 2.0
    tgt_df = hits_c.groupby("target").agg(
        mean_qmid=("qmid", "mean"),
        n_hits=("query", "count"),
        best_evalue=("evalue", "min"),
    ).reset_index()
    tgt_df["norm_mid"] = tgt_df["mean_qmid"] / qlen

    n_nterm = int((tgt_df["norm_mid"] < NTERM_CUTOFF).sum())
    n_cterm = int((tgt_df["norm_mid"] > CTERM_CUTOFF).sum())
    n_bridge = int(((tgt_df["norm_mid"] >= NTERM_CUTOFF) &
                     (tgt_df["norm_mid"] <= CTERM_CUTOFF)).sum())

    if n_nterm < MIN_DOMAIN_TARGETS or n_cterm < MIN_DOMAIN_TARGETS:
        return None
    if n_bridge >= min(n_nterm, n_cterm):
        return None

    # Classify each hit's target
    tgt_class = {}
    for _, row in tgt_df.iterrows():
        nm = row["norm_mid"]
        if nm < NTERM_CUTOFF:
            tgt_class[row["target"]] = 1
        elif nm > CTERM_CUTOFF:
            tgt_class[row["target"]] = 2
        else:
            tgt_class[row["target"]] = 0
    hits_c["domain"] = hits_c["target"].map(tgt_class).fillna(0).astype(int)

    nterm_hits = hits_c[hits_c["domain"] == 1]
    cterm_hits = hits_c[hits_c["domain"] == 2]
    nterm_ivs = merge_intervals(
        list(zip(nterm_hits["qstart"].tolist(), nterm_hits["qend"].tolist()))
    ) if len(nterm_hits) else []
    cterm_ivs = merge_intervals(
        list(zip(cterm_hits["qstart"].tolist(), cterm_hits["qend"].tolist()))
    ) if len(cterm_hits) else []

    nterm_cov = intervals_coverage(nterm_ivs, qlen)
    cterm_cov = intervals_coverage(cterm_ivs, qlen)
    total_cov = intervals_coverage(
        merge_intervals(nterm_ivs + cterm_ivs), qlen)

    nterm_range = (f"{nterm_ivs[0][0]}-{nterm_ivs[-1][1]}"
                   if nterm_ivs else "")
    cterm_range = (f"{cterm_ivs[0][0]}-{cterm_ivs[-1][1]}"
                   if cterm_ivs else "")

    # Best target subject coverage per domain
    def _best_target_tcov(domain_hits: pd.DataFrame) -> float:
        """Find the best target (highest total bits) and compute its tcov."""
        if domain_hits.empty:
            return 0.0
        tbits = domain_hits.groupby("target")["bits"].sum()
        best_tgt = tbits.idxmax()
        th = domain_hits[domain_hits["target"] == best_tgt]
        tlen = int(th["tlen"].iloc[0])
        tivs = merge_intervals(
            list(zip(th["tstart"].tolist(), th["tend"].tolist())))
        return intervals_coverage(tivs, tlen)

    nterm_best_tcov = _best_target_tcov(nterm_hits)
    cterm_best_tcov = _best_target_tcov(cterm_hits)

    return {
        "query_id": query_id, "qlen": qlen,
        "n_nterm_targets": n_nterm,
        "n_cterm_targets": n_cterm,
        "n_bridge_targets": n_bridge,
        "nterm_qcov": round(nterm_cov, 4),
        "cterm_qcov": round(cterm_cov, 4),
        "total_qcov": round(total_cov, 4),
        "nterm_best_tcov": round(nterm_best_tcov, 4),
        "cterm_best_tcov": round(cterm_best_tcov, 4),
        "nterm_range": nterm_range,
        "cterm_range": cterm_range,
    }, hits_c


# ---------------------------------------------------------------------------
# Stage 2: Streaming Pattern B scan
# ---------------------------------------------------------------------------

def streaming_pattern_b_scan(mmseqs_path: str, evalue: float,
                              min_hits: int, min_coverage: float,
                              summary_path: str, details_path: str,
                              chunk_size: int = 2_000_000,
                              verbose: bool = False) -> pd.DataFrame:
    """Scan full MMseqs2 output for Pattern B candidates.

    Reads the large MMseqs2 file in chunks with e-value filtering, accumulates
    all filtered hits, then analyses each query for bimodal target distribution.
    Applies min(nterm_best_tcov, cterm_best_tcov) > min_coverage filter.
    """
    print("[Pattern B streaming] Scanning full proteome ...", flush=True)
    t0 = time.time()

    if (os.path.isfile(summary_path) and os.path.getsize(summary_path) > 0):
        print(f"  Reusing existing result: {summary_path}")
        df = pd.read_csv(summary_path, sep="\t")
        print(f"  {len(df)} candidates (best target tcov > {min_coverage})")
        return df

    # --- Phase 1: chunked read with e-value filtering ---
    has_hdr = detect_header(mmseqs_path)
    reader = pd.read_csv(
        mmseqs_path, sep="\t",
        header=0 if has_hdr else None,
        names=None if has_hdr else MMSEQS_COLUMNS,
        usecols=COLS_HIT,
        dtype={c: MMSEQS_DTYPES[c] for c in COLS_HIT},
        chunksize=chunk_size,
    )

    parts: List[pd.DataFrame] = []
    n_chunks = 0
    n_rows_raw = 0
    for chunk in reader:
        n_chunks += 1
        n_rows_raw += len(chunk)
        filtered = chunk[chunk["evalue"] <= evalue]
        if len(filtered) > 0:
            parts.append(filtered)
        if verbose and n_chunks % 10 == 0:
            n_filt = sum(len(p) for p in parts)
            print(f"    chunk {n_chunks}: {n_rows_raw:,} rows read, "
                  f"{n_filt:,} pass e-value filter", flush=True)

    if not parts:
        print("  No hits pass e-value filter.")
        pd.DataFrame(columns=[
            "query_id", "qlen", "n_nterm_targets", "n_cterm_targets",
            "n_bridge_targets", "nterm_qcov", "cterm_qcov", "total_qcov",
            "nterm_best_tcov", "cterm_best_tcov",
            "nterm_range", "cterm_range",
        ]).to_csv(summary_path, sep="\t", index=False)
        return pd.DataFrame()

    all_hits = pd.concat(parts, ignore_index=True)
    del parts
    elapsed_read = time.time() - t0
    print(f"  Read {n_rows_raw:,} rows in {n_chunks} chunks ({elapsed_read:.1f}s), "
          f"{len(all_hits):,} pass e-value filter "
          f"({all_hits['query'].nunique():,} queries)")

    # --- Phase 2: per-query Pattern B analysis ---
    summaries: List[dict] = []
    detail_frames: List[pd.DataFrame] = []
    n_queries = 0

    for qid, grp in all_hits.groupby("query"):
        n_queries += 1
        result = _analyze_single_query_b(qid, grp, min_hits)
        if result is not None:
            summaries.append(result[0])
            detail_frames.append(result[1])

    elapsed = time.time() - t0
    n_before_cov = len(summaries)
    print(f"  Analysed {n_queries:,} queries ({elapsed:.1f}s)")
    print(f"  Found {n_before_cov} Pattern B candidates (before coverage filter)")

    # --- Phase 3: subject coverage filter ---
    if summaries:
        summary = pd.DataFrame(summaries)
        summary["min_tcov"] = summary[["nterm_best_tcov", "cterm_best_tcov"]].min(axis=1)
        summary = summary[summary["min_tcov"] > min_coverage]
        summary = summary.sort_values("min_tcov", ascending=False)
        summary = summary.drop(columns=["min_tcov"])
        print(f"  After min(nterm_best_tcov, cterm_best_tcov) > {min_coverage}: "
              f"{len(summary)} candidates")

        if len(summary) > 0:
            # Keep only detail frames for surviving candidates
            keep_ids = set(summary["query_id"])
            details = pd.concat(
                [d for d in detail_frames
                 if d["query"].iloc[0] in keep_ids],
                ignore_index=True)
            summary.to_csv(summary_path, sep="\t", index=False)
            details.to_csv(details_path, sep="\t", index=False)
            print(f"  {summary_path}")
            print(f"  {details_path}")
            return summary

    pd.DataFrame(columns=[
        "query_id", "qlen", "n_nterm_targets", "n_cterm_targets",
        "n_bridge_targets", "nterm_qcov", "cterm_qcov", "total_qcov",
        "nterm_range", "cterm_range",
    ]).to_csv(summary_path, sep="\t", index=False)
    return pd.DataFrame()


# ---------------------------------------------------------------------------
# Stage 3: Self-blast (candidates vs reference, self-excluded)
# ---------------------------------------------------------------------------

def run_self_blast_pattern_b(
    cand_fasta: str, ref_fasta: str,
    summary_path: str, details_path: str,
    mmseqs_path: str, threads: int, tmp_dir: str,
    evalue: float, min_hits: int, min_coverage: float,
    verbose: bool = False,
) -> pd.DataFrame:
    """Run candidates vs reference proteome (self-excluded), Pattern B."""
    print("[Self-blast] Candidates vs reference proteome ...", flush=True)
    t0 = time.time()

    if (os.path.isfile(summary_path) and os.path.getsize(summary_path) > 0):
        print(f"  Reusing existing result: {summary_path}")
        return pd.read_csv(summary_path, sep="\t")

    # MMseqs2 search
    raw_tsv = summary_path.replace(".summary.tsv", ".raw.tsv")
    os.makedirs(tmp_dir, exist_ok=True)
    run_mmseqs_search(cand_fasta, ref_fasta, raw_tsv,
                      mmseqs_path, threads, tmp_dir,
                      evalue=1e-3, verbose=verbose)

    # Read and filter
    hits = read_mmseqs_result(raw_tsv, evalue, exclude_self=True)
    print(f"  {len(hits):,} hits (self-excluded) for "
          f"{hits['query'].nunique()} candidates")

    # Pattern B analysis
    summaries: List[dict] = []
    detail_frames: List[pd.DataFrame] = []
    for qid, grp in hits.groupby("query"):
        result = _analyze_single_query_b(qid, grp, min_hits)
        if result is not None:
            summaries.append(result[0])
            detail_frames.append(result[1])

    elapsed = time.time() - t0
    print(f"  Found {len(summaries)} Pattern B candidates "
          f"(before coverage filter, {elapsed:.1f}s)")

    if summaries:
        summary = pd.DataFrame(summaries)
        summary["min_tcov"] = summary[["nterm_best_tcov", "cterm_best_tcov"]].min(axis=1)
        summary = summary[summary["min_tcov"] > min_coverage]
        summary = summary.sort_values("min_tcov", ascending=False)
        summary = summary.drop(columns=["min_tcov"])
        print(f"  After min(nterm_best_tcov, cterm_best_tcov) > {min_coverage}: "
              f"{len(summary)} candidates")

        if len(summary) > 0:
            keep_ids = set(summary["query_id"])
            details = pd.concat(
                [d for d in detail_frames
                 if d["query"].iloc[0] in keep_ids],
                ignore_index=True)
            summary.to_csv(summary_path, sep="\t", index=False)
            details.to_csv(details_path, sep="\t", index=False)
            print(f"  {summary_path}")
            print(f"  {details_path}")
            return summary

    pd.DataFrame(columns=[
        "query_id", "qlen", "n_nterm_targets", "n_cterm_targets",
        "n_bridge_targets", "nterm_qcov", "cterm_qcov", "total_qcov",
        "nterm_best_tcov", "cterm_best_tcov",
        "nterm_range", "cterm_range",
    ]).to_csv(summary_path, sep="\t", index=False)
    return pd.DataFrame()


# ---------------------------------------------------------------------------
# Stage 4: SwissProt description annotation (enhanced)
# ---------------------------------------------------------------------------

def _parse_sp_description(theader: str) -> str:
    """Extract protein name from UniProt FASTA header."""
    parts = str(theader).split(" ", 1)
    if len(parts) < 2:
        return str(theader)
    desc = parts[1]
    idx = desc.find(" OS=")
    if idx > 0:
        desc = desc[:idx]
    return desc.strip()


def _parse_sp_gene(theader: str) -> str:
    """Extract gene name from UniProt header GN= field."""
    m = re.search(r"GN=(\S+)", str(theader))
    return m.group(1) if m else ""


def _count_nonoverlap_sp_proteins(grp: pd.DataFrame,
                                   overlap_frac: float = 0.50
                                   ) -> List[dict]:
    """Greedy non-overlapping protein tiling of a query.

    For each UniProt target, merge its query alignment intervals.
    Then greedily select targets by bitscore, skipping any target whose
    query region overlaps >overlap_frac with already-selected targets.

    Returns list of selected proteins sorted by qstart:
      [{"target": ..., "gene": ..., "desc": ..., "qstart": ..., "qend": ...,
        "bits": ..., "tcov": ...}, ...]
    """
    # Per-target: merged query intervals, target intervals, best bits, header
    target_info: Dict[str, dict] = {}
    for target, tgrp in grp.groupby("target"):
        q_ivs = merge_intervals(
            list(zip(tgrp["qstart"].tolist(), tgrp["qend"].tolist())))
        t_ivs = merge_intervals(
            list(zip(tgrp["tstart"].tolist(), tgrp["tend"].tolist())))
        best_row = tgrp.loc[tgrp["bits"].idxmax()]
        tlen = int(best_row["tlen"])
        tcov = intervals_coverage(t_ivs, tlen)
        target_info[target] = {
            "intervals": q_ivs,
            "qstart": q_ivs[0][0],
            "qend": q_ivs[-1][1],
            "span": sum(e - s for s, e in q_ivs),
            "bits": float(best_row["bits"]),
            "theader": best_row["theader"],
            "tcov": tcov,
        }

    # Sort by bitscore descending
    sorted_targets = sorted(target_info, key=lambda t: target_info[t]["bits"],
                            reverse=True)

    selected: List[dict] = []
    covered: List[Tuple[int, int]] = []

    for target in sorted_targets:
        info = target_info[target]
        # Compute overlap with already-covered regions
        ov = 0
        for s, e in info["intervals"]:
            for cs, ce in covered:
                os_, oe = max(s, cs), min(e, ce)
                if oe > os_:
                    ov += oe - os_
        if info["span"] > 0 and ov / info["span"] > overlap_frac:
            continue  # too much overlap, skip

        selected.append({
            "target": target,
            "gene": _parse_sp_gene(info["theader"]),
            "desc": _parse_sp_description(info["theader"]),
            "qstart": info["qstart"],
            "qend": info["qend"],
            "bits": info["bits"],
            "tcov": round(info["tcov"], 4),
        })
        covered = merge_intervals(covered + info["intervals"])

    # Sort by query position
    selected.sort(key=lambda x: x["qstart"])
    return selected


def annotate_swissprot_descriptions(
    cand_fasta: str, sp_fasta: str, out_tsv: str,
    mmseqs_path: str, threads: int, tmp_dir: str,
    evalue: float, verbose: bool = False,
) -> Dict[str, dict]:
    """Search candidates vs SwissProt FASTA, return per-query annotation.

    swissprot_n_genes = number of non-overlapping best-match UniProt
    proteins that tile the query (greedy by bitscore).
    Also reports N-term/C-term best hits with query start/end.
    """
    print("[SwissProt annotation] Searching vs UniProt SwissProt ...")
    os.makedirs(tmp_dir, exist_ok=True)

    run_mmseqs_search(cand_fasta, sp_fasta, out_tsv,
                      mmseqs_path, threads, tmp_dir,
                      evalue=1e-3, verbose=verbose,
                      fmt=",".join(SWISSPROT_COLS))

    df = pd.read_csv(out_tsv, sep="\t", header=None, names=SWISSPROT_COLS)
    df = df[df["evalue"] <= evalue]

    results: Dict[str, dict] = {}
    for qid, grp in df.groupby("query"):
        qlen = int(grp["qlen"].iloc[0])

        # --- Non-overlapping protein tiling ---
        selected = _count_nonoverlap_sp_proteins(grp)
        n_genes = len(selected)

        # Best overall hit
        best = grp.loc[grp["bits"].idxmax()]
        entry: dict = {
            "swissprot_n_genes": n_genes,
            "best_sp_hit": best["target"],
            "best_sp_gene": _parse_sp_gene(best["theader"]),
            "best_sp_evalue": best["evalue"],
        }

        # Min subject coverage across non-overlapping proteins
        if selected:
            entry["sp_tcov"] = round(min(s["tcov"] for s in selected), 4)
        else:
            entry["sp_tcov"] = 0.0

        # Detail string: gene(qstart-qend,tcov) for each non-overlapping protein
        detail_parts = []
        for sel in selected:
            g = sel["gene"] or sel["target"]
            detail_parts.append(f"{g}({sel['qstart']}-{sel['qend']},tcov={sel['tcov']:.2f})")
        entry["sp_genes_detail"] = "; ".join(detail_parts)

        # N-term / C-term best hits (based on query midpoint classification)
        grp_c = grp.copy()
        grp_c["qmid_norm"] = ((grp_c["qstart"] + grp_c["qend"]) / 2.0) / qlen
        nterm = grp_c[grp_c["qmid_norm"] < NTERM_CUTOFF]
        cterm = grp_c[grp_c["qmid_norm"] > CTERM_CUTOFF]

        if len(nterm) > 0:
            nbest = nterm.loc[nterm["bits"].idxmax()]
            entry["nterm_sp_hit"] = nbest["target"]
            entry["nterm_sp_gene"] = _parse_sp_gene(nbest["theader"])
            entry["nterm_qstart"] = int(nterm["qstart"].min())
            entry["nterm_qend"] = int(nterm["qend"].max())
            entry["nterm_sp_desc"] = _parse_sp_description(nbest["theader"])
        else:
            entry["nterm_sp_hit"] = ""
            entry["nterm_sp_gene"] = ""
            entry["nterm_qstart"] = ""
            entry["nterm_qend"] = ""
            entry["nterm_sp_desc"] = ""

        if len(cterm) > 0:
            cbest = cterm.loc[cterm["bits"].idxmax()]
            entry["cterm_sp_hit"] = cbest["target"]
            entry["cterm_sp_gene"] = _parse_sp_gene(cbest["theader"])
            entry["cterm_qstart"] = int(cterm["qstart"].min())
            entry["cterm_qend"] = int(cterm["qend"].max())
            entry["cterm_sp_desc"] = _parse_sp_description(cbest["theader"])
        else:
            entry["cterm_sp_hit"] = ""
            entry["cterm_sp_gene"] = ""
            entry["cterm_qstart"] = ""
            entry["cterm_qend"] = ""
            entry["cterm_sp_desc"] = ""

        entry["best_sp_desc"] = _parse_sp_description(best["theader"])

        results[qid] = entry

    print(f"  Annotated {len(results)} candidates with SwissProt descriptions")
    n_multi = sum(1 for v in results.values() if v["swissprot_n_genes"] >= 2)
    print(f"  swissprot_n_genes >= 2: {n_multi} candidates")
    return results


# ---------------------------------------------------------------------------
# Combined summary
# ---------------------------------------------------------------------------

def build_combined_summary(
    candidates: Set[str],
    pattern_b_sp_summary: pd.DataFrame,
    pattern_b_self_summary: pd.DataFrame,
    sp_annot: Dict[str, dict],
    output_path: str,
) -> pd.DataFrame:
    """Build combined summary: Pattern B (SwissProt + self-blast) + SP annotation.

    Column order: identifiers / evidence -> counts -> coordinates -> descriptions.
    """
    # Index Pattern B summaries for lookup
    b_sp_dict: Dict[str, dict] = {}
    if pattern_b_sp_summary is not None and len(pattern_b_sp_summary) > 0:
        for _, r in pattern_b_sp_summary.iterrows():
            b_sp_dict[r["query_id"]] = r.to_dict()

    b_self_dict: Dict[str, dict] = {}
    if pattern_b_self_summary is not None and len(pattern_b_self_summary) > 0:
        for _, r in pattern_b_self_summary.iterrows():
            b_self_dict[r["query_id"]] = r.to_dict()

    rows = []
    for cand in sorted(candidates):
        # Skip candidates without 2+ distinct non-overlapping proteins
        sa = sp_annot.get(cand, {}) if sp_annot else {}
        n_genes = sa.get("swissprot_n_genes", 0)
        if n_genes < 2:
            continue

        row: dict = {"candidate": cand}

        # Pattern B info (from Stage 2 — vs SwissProt)
        pb = b_sp_dict.get(cand, {})
        row["qlen"] = pb.get("qlen", "")
        row["nterm_best_tcov"] = pb.get("nterm_best_tcov", "")
        row["cterm_best_tcov"] = pb.get("cterm_best_tcov", "")

        # SwissProt annotation
        row["swissprot_n_genes"] = n_genes
        row["sp_tcov"] = sa.get("sp_tcov", "")

        # Evidence flags
        row["pattern_b_vs_swissprot"] = cand in b_sp_dict
        row["pattern_b_vs_self"] = cand in b_self_dict

        # Best overall hit
        row["best_sp_hit"] = sa.get("best_sp_hit", "")
        row["best_sp_gene"] = sa.get("best_sp_gene", "")
        row["best_sp_evalue"] = sa.get("best_sp_evalue", "")

        # N-term and C-term coordinates together
        row["nterm_sp_hit"] = sa.get("nterm_sp_hit", "")
        row["nterm_sp_gene"] = sa.get("nterm_sp_gene", "")
        row["nterm_qstart"] = sa.get("nterm_qstart", "")
        row["nterm_qend"] = sa.get("nterm_qend", "")
        row["cterm_sp_hit"] = sa.get("cterm_sp_hit", "")
        row["cterm_sp_gene"] = sa.get("cterm_sp_gene", "")
        row["cterm_qstart"] = sa.get("cterm_qstart", "")
        row["cterm_qend"] = sa.get("cterm_qend", "")

        # All non-overlapping proteins detail
        row["sp_genes_detail"] = sa.get("sp_genes_detail", "")

        # Descriptions last
        row["best_sp_desc"] = sa.get("best_sp_desc", "")
        row["nterm_sp_desc"] = sa.get("nterm_sp_desc", "")
        row["cterm_sp_desc"] = sa.get("cterm_sp_desc", "")

        rows.append(row)

    if not rows:
        print("  No candidates with swissprot_n_genes >= 2")
        pd.DataFrame().to_csv(output_path, sep="\t", index=False)
        return pd.DataFrame()

    summary = pd.DataFrame(rows)

    evidence_cols = ["pattern_b_vs_swissprot", "pattern_b_vs_self"]
    summary["evidence_count"] = (
        summary[evidence_cols].sum(axis=1).astype(int))

    # Sort: evidence_count desc, then swissprot_n_genes desc
    summary = summary.sort_values(
        ["evidence_count", "swissprot_n_genes"], ascending=[False, False])

    # Final column order: ID -> tcov -> evidence -> hits/coords -> detail -> descriptions
    col_order = [
        "candidate", "qlen", "nterm_best_tcov", "cterm_best_tcov", "sp_tcov",
        "swissprot_n_genes",
        "pattern_b_vs_swissprot", "pattern_b_vs_self", "evidence_count",
        "best_sp_hit", "best_sp_gene", "best_sp_evalue",
        "nterm_sp_hit", "nterm_sp_gene",
        "nterm_qstart", "nterm_qend",
        "cterm_sp_hit", "cterm_sp_gene",
        "cterm_qstart", "cterm_qend",
        "sp_genes_detail",
        "best_sp_desc", "nterm_sp_desc", "cterm_sp_desc",
    ]
    summary = summary[[c for c in col_order if c in summary.columns]]

    summary.to_csv(output_path, sep="\t", index=False)
    n_sources = len(evidence_cols)
    print(f"  Evidence sources ({n_sources}): " +
          ", ".join(evidence_cols))
    print(f"  Total candidates (swissprot_n_genes >= 2): {len(summary)}")
    for threshold in range(n_sources, 0, -1):
        n = (summary["evidence_count"] >= threshold).sum()
        print(f"  >= {threshold}/{n_sources} evidence: {n} candidates")
    if "swissprot_n_genes" in summary.columns:
        n2 = (summary["swissprot_n_genes"] >= 2).sum()
        n3 = (summary["swissprot_n_genes"] > 2).sum()
        print(f"  swissprot_n_genes >= 2: {n2} candidates")
        print(f"  swissprot_n_genes >  2: {n3} candidates")
    print(f"  {output_path}")
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Standalone whole-proteome chimera scan: "
                    "MMseqs2 search -> Pattern B -> self-blast -> "
                    "SwissProt annotation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Required
    p.add_argument("--ref-fasta", required=True,
                   help="Reference proteome FASTA")
    p.add_argument("--swissprot-db", required=True,
                   help="MMseqs2 SwissProt DB or FASTA (search target)")

    # Optional inputs
    p.add_argument("--swissprot-fasta",
                   help="UniProt SwissProt FASTA with descriptions "
                        "(for Stage 4; skipped if not provided)")
    p.add_argument("--mmseqs-out",
                   help="Pre-existing MMseqs2 result TSV "
                        "(skip Stage 1 if provided)")

    # MMseqs2 settings
    p.add_argument("--mmseqs-path", default="mmseqs",
                   help="Path to mmseqs binary")
    p.add_argument("--mmseqs-threads", type=int, default=16)

    # Analysis params
    p.add_argument("--evalue", type=float, default=1e-5,
                   help="E-value cutoff for Pattern B filtering")
    p.add_argument("--min-hits", type=int, default=5,
                   help="Min hits per query for Pattern B")
    p.add_argument("--min-coverage", type=float, default=0.80,
                   help="Min best-target subject coverage per domain")
    p.add_argument("--chunk-size", type=int, default=2_000_000,
                   help="Rows per chunk for streaming scan")

    # Output
    p.add_argument("--output-prefix", default="chimera_fullscan")
    p.add_argument("--output-dir", default=".")
    p.add_argument("--verbose", action="store_true")

    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)

    prefix = os.path.join(args.output_dir, args.output_prefix)

    # ==================================================================
    # Stage 1: MMseqs2 search (reference vs SwissProt)
    # ==================================================================
    print("=" * 60)
    print("STAGE 1: MMseqs2 search -- reference vs SwissProt")
    print("=" * 60)

    mmseqs_out = args.mmseqs_out
    if mmseqs_out is None:
        mmseqs_out = f"{prefix}.mmseqs.tsv"
        tmp_dir = os.path.join(args.output_dir, "mmseqs_tmp_fullscan")
        os.makedirs(tmp_dir, exist_ok=True)
        run_mmseqs_search(
            args.ref_fasta, args.swissprot_db, mmseqs_out,
            args.mmseqs_path, args.mmseqs_threads, tmp_dir,
            evalue=1e-3, verbose=args.verbose)
    else:
        print(f"  Using pre-existing MMseqs2 result: {mmseqs_out}")

    file_size_gb = os.path.getsize(mmseqs_out) / (1024**3)
    print(f"  Result file: {mmseqs_out} ({file_size_gb:.2f} GB)\n")

    # ==================================================================
    # Stage 2: Streaming Pattern B scan (coverage > min_coverage)
    # ==================================================================
    print("=" * 60)
    print(f"STAGE 2: Pattern B scan (best target tcov > {args.min_coverage})")
    print("=" * 60)

    summary_path = f"{prefix}.pattern_b.summary.tsv"
    details_path = f"{prefix}.pattern_b.details.tsv"

    pattern_b_summary = streaming_pattern_b_scan(
        mmseqs_out, args.evalue, args.min_hits, args.min_coverage,
        summary_path, details_path,
        chunk_size=args.chunk_size, verbose=args.verbose)

    if pattern_b_summary is None or len(pattern_b_summary) == 0:
        print("  No Pattern B candidates found. Pipeline complete.\n")
        return

    candidates = set(pattern_b_summary["query_id"])
    print(f"  {len(candidates)} candidates for downstream analysis\n")

    # ==================================================================
    # Stage 3: Self-blast (candidates vs reference, self-excluded)
    # ==================================================================
    print("=" * 60)
    print("STAGE 3: Self-blast -- candidates vs reference (self-excluded)")
    print("=" * 60)

    # Extract candidate sequences
    cand_fasta = f"{prefix}.candidates.fasta"
    n = extract_fasta_subset(args.ref_fasta, candidates, cand_fasta)
    print(f"  Extracted {n} candidate sequences -> {cand_fasta}")

    self_summary_path = f"{prefix}.self_blast.summary.tsv"
    self_details_path = f"{prefix}.self_blast.details.tsv"
    self_tmp = os.path.join(args.output_dir, "mmseqs_tmp_self_blast")

    pattern_b_self = run_self_blast_pattern_b(
        cand_fasta, args.ref_fasta,
        self_summary_path, self_details_path,
        args.mmseqs_path, args.mmseqs_threads, self_tmp,
        args.evalue, args.min_hits, args.min_coverage,
        verbose=args.verbose)
    print()

    # ==================================================================
    # Stage 4: SwissProt description annotation (enhanced)
    # ==================================================================
    sp_annot: Dict[str, dict] = {}
    if args.swissprot_fasta:
        print("=" * 60)
        print("STAGE 4: SwissProt description annotation")
        print("=" * 60)

        sp_out = f"{prefix}.vs_swissprot_described.tsv"
        sp_tmp = os.path.join(args.output_dir, "mmseqs_tmp_sp_desc_fs")
        sp_annot = annotate_swissprot_descriptions(
            cand_fasta, args.swissprot_fasta, sp_out,
            args.mmseqs_path, args.mmseqs_threads, sp_tmp,
            args.evalue, verbose=args.verbose)
        print()
    else:
        print("\n  --swissprot-fasta not provided, skipping Stage 4.\n")

    # ==================================================================
    # Combined summary
    # ==================================================================
    print("=" * 60)
    print("COMBINED SUMMARY")
    print("=" * 60)

    combined_path = f"{prefix}.combined_summary.tsv"
    build_combined_summary(
        candidates, pattern_b_summary, pattern_b_self,
        sp_annot, combined_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
