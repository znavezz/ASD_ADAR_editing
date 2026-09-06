#!/usr/bin/env python
"""Build the committed test fixture: a handful of real variants plus their real VEP output.

The pipeline's inputs are a 463 MB VCF, a 3.2 GB genome and a 144 MB VEP output file, none of
which can be committed or used in CI. This cuts a tiny, deterministic slice that still
exercises every shape the code branches on, so Phases 1-2 can be refactored against real data
without the network, Docker, or 12 GB of resources.

Selection is by *criterion*, not by hand-picked accession: each case is defined as a predicate
and the lowest-id matching variant is taken, so the set is reproducible and its rationale is
self-documenting. Every chosen variant must appear in the real VEP output, and the matching
VEP lines are sliced out verbatim -- including the quirks of the `Extra` column, which is the
whole point.

Usage::

    python scripts/build_fixtures.py            # write tests/data/fixtures/
    python scripts/build_fixtures.py --check    # fail if the committed fixture is stale
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SNAPSHOT = REPO_ROOT / "Output" / "Intermediate_tables" / "00_db_snapshot" / "df_cds.csv"
VEP_OUTPUT = REPO_ROOT / "Output" / "VEP" / "varicarta_vepped.txt"
TARGET_DIR = REPO_ROOT / "tests" / "data" / "fixtures"

MAX_BYTES = 1_000_000  # fixtures must stay tiny; CI fails above 5x this

#: Variants that cannot come from the CDS snapshot because it is SNV-only. Taken as the
#: lowest-id Deletion and Insertion in the database, so they are still deterministic.
NON_SNV = [
    ("1:103616:CTTTTTTTT:CTTTT", "deletion", "class != SNV; exercises the CADD SNV-only gate"),
    ("1:840685:G:GT", "insertion", "class != SNV; ref/alt lengths differ, end != start"),
]


def criteria(df: pd.DataFrame) -> dict[str, tuple[pd.Series, str]]:
    """Each entry: name -> (row mask, why this case matters)."""
    plus, minus = df.gene_strand == "+", df.gene_strand == "-"
    # NOTE: REF/ALT in the snapshot are strand-CORRECTED (coding strand);
    # REF.VEP/ALT.VEP are the raw genomic alleles. See docs/provenance.md.
    genomic_g, genomic_a = df["REF.VEP"] == "G", df["ALT.VEP"] == "A"
    genomic_c, genomic_t = df["REF.VEP"] == "C", df["ALT.VEP"] == "T"
    alt_has_a = df[["alt_codon_nt1", "alt_codon_nt2", "alt_codon_nt3"]].eq("A").any(axis=1)
    clin = df.CLIN_SIG.fillna("")

    return {
        "g2a_plus_strand": (
            genomic_g & genomic_a & plus,
            "directly correctable, plus strand: genomic G>A == coding G>A",
        ),
        "c2t_minus_strand": (
            genomic_c & genomic_t & minus,
            "directly correctable, minus strand: genomic C>T == coding G>A. "
            "The case the 218-vs-138 error got wrong",
        ),
        "g2a_on_minus_not_correctable": (
            genomic_g & genomic_a & minus,
            "NOT correctable despite being genomic G>A: on a minus-strand gene it reads "
            "as C>T. The complement of the case above; both are needed to pin the strand rule",
        ),
        "stopgained": (
            df.New_Consequence == "Nonsense",
            "the Rescue path: stop codon recoded to TGG",
        ),
        "missense_alt_codon_has_A": (
            (df.New_Consequence == "Missense") & alt_has_a,
            "the Improve path: an editable adenosine inside the mutant codon",
        ),
        # The criterion above is satisfied by a variant that is *also* directly correctable,
        # so deduplication folds the two together and the Improve path ends up covered by a
        # Fix case. This one excludes revertible variants, so the fixture contains a variant
        # that residue optimisation actually classifies.
        "improve_not_also_fix": (
            (df.New_Consequence == "Missense")
            & alt_has_a
            & ~((genomic_g & genomic_a & plus) | (genomic_c & genomic_t & minus)),
            "residue optimisation with no simpler route: missense, editable adenosine in the "
            "mutant codon, and NOT directly revertible",
        ),
        "splice": (df.New_Consequence == "Splice", "splice-site consequence bucket"),
        "synonymous": (df.New_Consequence == "Synonymous", "in the CDS set but not treatable"),
        "start_lost": (df.New_Consequence == "Start_lost", "rare consequence class (n=10)"),
        "multi_consequence": (
            df.Consequence.str.contains(",", na=False),
            "several VEP terms on one variant; exercises consequence collapsing and "
            "most-severe ranking",
        ),
        "missing_sift": (
            df.SIFT_score.isna() & (df.New_Consequence == "Missense"),
            "missense with no SIFT score; must degrade, not crash",
        ),
        "chrX": (df.chr.astype(str) == "X", "non-autosome; guards chromosome handling"),
        "clinvar_pathogenic": (
            clin.str.contains("pathogenic", case=False)
            & ~clin.str.contains("conflict|benign", case=False),
            "ClinVar pathogenic under the Figure 1C rule",
        ),
        "codon_position_1": (df.position_in_codon == 1, "codon offset boundary"),
        "codon_position_3": (df.position_in_codon == 3, "codon offset boundary"),
        "no_rs_id": (df.rs_id.isna(), "no Existing_variation; rs_id must be NULL not empty"),
    }


def select(df: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Return [(variant_id, criterion, rationale)], lowest id per criterion, deduplicated."""
    df = df.sort_values("id").reset_index(drop=True)
    vid = (
        df.chr.astype(str) + ":" + df.start.astype(str) + ":" + df["REF.VEP"] + ":" + df["ALT.VEP"]
    )
    chosen: dict[str, tuple[str, str]] = {}
    for name, (mask, why) in criteria(df).items():
        hits = vid[mask]
        if hits.empty:
            continue
        picked = hits.iloc[0]
        if picked in chosen:  # already covered; record the extra criterion it satisfies
            prev_crit, prev_why = chosen[picked]
            chosen[picked] = (f"{prev_crit} + {name}", f"{prev_why}; also: {why}")
        else:
            chosen[picked] = (name, why)
    out = [(v, c, w) for v, (c, w) in chosen.items()]
    out.extend((v, c, w) for v, c, w in NON_SNV)
    return out


def slice_vep(variant_ids: set[str]) -> tuple[list[str], list[str], set[str]]:
    """Stream the VEP output, keeping the ## header and rows for the chosen variants."""
    header, rows, found = [], [], set()
    with VEP_OUTPUT.open() as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            vid = line.split("\t", 1)[0]
            if vid in variant_ids:
                rows.append(line)
                found.add(vid)
    return header, rows, found


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="Fail if the committed fixture is stale.")
    args = ap.parse_args(argv)

    for path in (SNAPSHOT, VEP_OUTPUT):
        if not path.exists():
            raise SystemExit(f"Required input missing: {path}")

    picks = select(pd.read_csv(SNAPSHOT, low_memory=False))
    ids = {v for v, _, _ in picks}
    print(f"Selected {len(picks)} variants covering {len(criteria(pd.read_csv(SNAPSHOT, nrows=1)))} criteria")

    header, rows, found = slice_vep(ids)
    missing = ids - found
    if missing:
        raise SystemExit(
            "These variants are not in the VEP output, so the fixture would be "
            f"inconsistent: {sorted(missing)}"
        )

    picks.sort(key=lambda t: t[0])
    vcf = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    readme = [
        "# Test fixtures",
        "",
        "Real variants and their **real** VEP output, sliced from the paper run.",
        "Regenerate with `python scripts/build_fixtures.py`; verify with `--check`.",
        "",
        "Selection is by criterion (lowest matching variant id), never hand-picked, so this set",
        "is reproducible and each entry documents why it is here.",
        "",
        "`vep_slice.txt` keeps the original `##` header, so the VEP version and the bundled",
        "ClinVar/SIFT/dbSNP releases travel with the fixture (see `docs/provenance.md`).",
        "",
        "| variant (chr:pos:ref:alt) | criterion | why it is in the fixture |",
        "|---|---|---|",
    ]
    for v, crit, why in picks:
        chrom, pos, ref, alt = v.split(":")
        vcf.append(f"{chrom}\t{pos}\t{v}\t{ref}\t{alt}\t.\t.\t.")
        readme.append(f"| `{v}` | {crit} | {why} |")

    TARGET_DIR.mkdir(parents=True, exist_ok=True)
    artifacts = {
        "variants.vcf": "\n".join(vcf) + "\n",
        "vep_slice.txt": "".join(header) + "".join(sorted(rows)),
        "README.md": "\n".join(readme) + "\n",
    }

    if args.check:
        for name, content in artifacts.items():
            path = TARGET_DIR / name
            if not path.exists():
                raise SystemExit(f"{path.relative_to(REPO_ROOT)} missing; run without --check.")
            if path.read_text() != content:
                raise SystemExit(f"{path.relative_to(REPO_ROOT)} is stale; regenerate it.")
        print("Fixtures are current.")
        return 0

    total = 0
    for name, content in artifacts.items():
        (TARGET_DIR / name).write_text(content)
        total += len(content)
        print(f"  wrote {name:<16} {len(content):>8,} bytes")
    if total > MAX_BYTES:
        raise SystemExit(f"Fixture total {total:,} bytes exceeds the {MAX_BYTES:,} budget.")
    print(f"Total {total:,} bytes (budget {MAX_BYTES:,}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
