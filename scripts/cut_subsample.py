#!/usr/bin/env python
"""Cut a small, deterministic, coverage-forcing subsample of the CDS snapshot.

``Output/Intermediate_tables/00_db_snapshot/df_cds.csv`` is the 9,962-row table every
figure downstream of Figure 1 is built from, and it is gitignored. Committing a checksum
(see scripts/verify_counts.py) catches *whether* it changed but says nothing about *what*
changed. This cuts a ~200-row subsample, committed to the repo, so a refactor that alters
individual values -- not just counts -- fails a test with a readable diff.

The sample is chosen so every interesting shape is present: each consequence class, both
gene strands, the ADAR "Fix" patterns (G>A on +, C>T on -), rows missing SIFT, rows with
several consequence terms, NMD-escaping and not, and each SFARI gene score. The remainder
is a fixed stride over sorted ids, so the output is byte-stable across runs.

Usage::

    python scripts/cut_subsample.py            # writes tests/data/golden/cds_subsample.csv
    python scripts/cut_subsample.py --check    # fails if the committed file is stale
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE = REPO_ROOT / "Output" / "Intermediate_tables" / "00_db_snapshot" / "df_cds.csv"
TARGET = REPO_ROOT / "tests" / "data" / "golden" / "cds_subsample.csv"
TARGET_ROWS = 200


def _strata(df: pd.DataFrame) -> dict[str, pd.Series]:
    """Named boolean masks; the lowest id matching each is force-included."""
    fix_plus = (df["REF"] == "G") & (df["ALT"] == "A") & (df["gene_strand"] == "+")
    fix_minus = (df["REF"] == "C") & (df["ALT"] == "T") & (df["gene_strand"] == "-")
    strata: dict[str, pd.Series] = {
        "fix_pattern_plus_strand": fix_plus,
        "fix_pattern_minus_strand": fix_minus,
        "non_fix_pattern": ~(fix_plus | fix_minus),
        "multi_consequence": df["Consequence"].str.contains(",", na=False),
        "missing_sift": df["SIFT_score"].isna(),
        "has_sift": df["SIFT_score"].notna(),
        "missing_cds_position": df["cds_position"].isna(),
        "has_intron": df["intron"].notna(),
        "has_rs_id": df["rs_id"].notna(),
        "missing_rs_id": df["rs_id"].isna(),
        "has_clin_sig": df["CLIN_SIG"].notna(),
        "nmd_escaping": df["nmd_escaping_variant"].astype(bool),
        "not_nmd_escaping": ~df["nmd_escaping_variant"].astype(bool),
        "stop_codon_ref": df["ref_codon_stop"].astype(str).str.lower().eq("true"),
        "stop_codon_alt": df["alt_codon_stop"].astype(str).str.lower().eq("true"),
    }
    for value in sorted(df["New_Consequence"].dropna().unique()):
        strata[f"consequence={value}"] = df["New_Consequence"] == value
    for value in sorted(df["position_in_codon"].dropna().unique()):
        strata[f"position_in_codon={int(value)}"] = df["position_in_codon"] == value
    for value in sorted(df["sfari_gene_score"].dropna().unique()):
        strata[f"sfari_gene_score={int(value)}"] = df["sfari_gene_score"] == value
    for value in sorted(df["gene_strand"].dropna().unique()):
        strata[f"gene_strand={value}"] = df["gene_strand"] == value
    return strata


def build(df: pd.DataFrame, target_rows: int = TARGET_ROWS) -> pd.DataFrame:
    df = df.sort_values("id").reset_index(drop=True)
    keep: list[int] = []
    for mask in _strata(df).values():
        matching = df.loc[mask, "id"]
        if not matching.empty:
            keep.append(int(matching.iloc[0]))
    keep = sorted(set(keep))

    # Fill the remainder with an even stride over the sorted ids, so the sample spans the
    # whole table rather than clustering at whichever ids the strata happened to pick.
    remaining = target_rows - len(keep)
    if remaining > 0:
        pool = [i for i in df["id"].tolist() if i not in set(keep)]
        if pool:
            stride = max(1, len(pool) // remaining)
            keep.extend(pool[::stride][:remaining])
    return df[df["id"].isin(set(keep))].sort_values("id").reset_index(drop=True)


def coverage_report(sample: pd.DataFrame, full: pd.DataFrame) -> list[str]:
    """Strata present in the full table but absent from the sample."""
    sample_strata = _strata(sample)
    return [
        name
        for name, mask in _strata(full).items()
        if mask.any() and not sample_strata.get(name, pd.Series(dtype=bool)).any()
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true", help="Fail if the committed file is stale.")
    parser.add_argument("--rows", type=int, default=TARGET_ROWS)
    args = parser.parse_args(argv)

    if not SOURCE.exists():
        raise SystemExit(f"Source snapshot not found: {SOURCE}")

    full = pd.read_csv(SOURCE, low_memory=False)
    sample = build(full, args.rows)
    rendered = sample.to_csv(index=False)

    missing = coverage_report(sample, full)
    if missing:
        raise SystemExit("Sample does not cover: " + ", ".join(missing))

    if args.check:
        if not TARGET.exists():
            raise SystemExit(f"{TARGET.relative_to(REPO_ROOT)} is missing; run without --check.")
        if TARGET.read_text() != rendered:
            raise SystemExit(f"{TARGET.relative_to(REPO_ROOT)} is stale; regenerate it.")
        print(f"{TARGET.relative_to(REPO_ROOT)} is current ({len(sample)} rows).")
        return 0

    TARGET.parent.mkdir(parents=True, exist_ok=True)
    TARGET.write_text(rendered)
    print(
        f"Wrote {TARGET.relative_to(REPO_ROOT)}: {len(sample)} of {len(full)} rows, "
        f"{len(sample.columns)} columns, {len(rendered) / 1024:.0f} KB"
    )
    print(f"Covers {len(_strata(full))} strata.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
