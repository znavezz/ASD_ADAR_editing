#!/usr/bin/env python3
"""
Export variant supplementary tables (flat) from Hasura.

Two variant SETS are supported (select with --set; default both):
  * treatable — ADAR A>G-correctable SNVs in SFARI protein-coding genes with a
                StopGained/Missense/Splice consequence (the original query).
  * rescue    — the 1092 NonG2A StopGained SNVs: same SFARI/protein-coding filters
                but NOT directly G>A/C>T-correctable, consequence StopGained only.
                Rescue rows carry extra columns: rescue_sift and the resulting
                codon/aa (see below).

Each set is written in three formats (select with --format; default all):
  1. <set>.csv          — MERGED: one row per variant; all bystanders packed into
                          one `bystanders(...)` column as ` | `-joined tuples.
  2. <set>_wide.csv     — WIDE: one row per variant; each bystander in its own
                          column (bystander_1 … bystander_20, the +-10 nt window
                          max); unused slots are `NULL`.
  3. Supplementary_Tables[_rescue].xlsx — single-sheet workbook mirroring the merged
                          CSV (one row per variant; all bystanders packed into one
                          column). No delimiters -> renders cleanly in Google Sheets/Excel.

All three formats share the same variant column order. Positions are 1-based
(verified against hg19); position columns carry the base, e.g. `pos (1-based)`,
`cds_position (1-based)`.

Rescue resulting codon: the pipeline rescues a premature stop by ADAR read-through
to Trp — `_build_rescue_sift_query` in pipeline/phase3_guides.py always targets TGG, and
`rescue_sift` is the VEP SIFT of (original AA -> Trp). The resulting codon is NOT
stored in the DB but is constant by construction, so it is emitted here as
rescue_resulting_codon=TGG / rescue_resulting_aa=Trp (reproducible, no DB change).

Bystander tuple field order (merged + wide):
    consequence; genomic_position (1-based); sift_score; sift_prediction;
    cadd_phred; ref_codon; ref_aa; alt_codon; alt_aa
`consequence` is the single most-severe region for the bystander (CONSEQUENCE_RANK).

Query notes: SIFT + CADD_phred are fetched by name and split into columns; results
page via limit/offset; variants_guides is restricted to the `pre_mRNA` guide
(off-targets, rescue_sift and bystanders all come from it); bystander id +
genomic_position are fetched for the tidy detail sheet.

Usage:
    python scripts/export_tables.py [--set all|treatable|rescue]
        [--format all|merged|wide|xlsx] [--out-dir DIR] [--limit N] [--page-size N]
"""

import os
import csv
import time
import argparse
import logging
from pathlib import Path

import requests
from dotenv import load_dotenv
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

# ── Setup ──
PROJECT_ROOT = Path(__file__).resolve().parent

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

env_path = PROJECT_ROOT / ".env"
if not env_path.exists():
    raise RuntimeError("Could not find .env file")
load_dotenv(env_path, override=True)

# Rescue read-through target (see phase3_guides._build_rescue_sift_query, target "TGG").
RESCUE_RESULTING_CODON = "TGG"
RESCUE_RESULTING_AA = "Trp"


def gql(query, variables):
    """Minimal Hasura client mirroring helpers.gql (admin secret + retries)."""
    url = str(os.environ["HASURA_URL"])
    headers = {"x-hasura-admin-secret": str(os.environ["HASURA_ADMIN_SECRET"])}
    for attempt in range(5):
        r = requests.post(url, json={"query": query, "variables": variables},
                          headers=headers, timeout=180)
        if r.status_code == 200:
            payload = r.json()
            if "errors" not in payload:
                return payload["data"]
            raise RuntimeError(payload["errors"])
        if r.status_code in (502, 503, 504):
            time.sleep(0.5 * (2 ** attempt))
            continue
        raise RuntimeError(f"HTTP {r.status_code}: {r.text[:800]}")
    raise RuntimeError("Hasura request failed after retries")


# ── Variant set WHERE clauses (only difference between the two sets) ──
WHERE_TREATABLE = """
      class: {_eq: "SNV"}
      _or: [
        {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
        {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
      ]
      variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}}
      gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}}
      variants_consequences: {consequence: {name: {_in: ["StopGained", "Missense", "Splice"]}}}
"""

WHERE_RESCUE = """
      class: {_eq: "SNV"},
      nmd_escaping_variant: {_eq: false},
      _not: {
        _or: [
          {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
          {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
        ]
      }
      variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}}
      gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}}
      variants_consequences: {consequence: {name: {_eq: "StopGained"}}}

"""

# Missense-improvement (neighbor-drill "Improve") set: the 158 reportable variants —
# a deleterious residue (variant SIFT < 0.05) whose best neighbor A->G edit is tolerated
# (SIFT >= 0.05) & better. Extra columns describe the chosen edit (best neighbor edit).
WHERE_IMPROVE = """
      class: {_eq: "SNV"}
      _not: {
        _or: [
          {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
          {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
        ]
      }
      variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}}
      gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}}
      variants_consequences: {consequence: {name: {_eq: "Missense"}}}
      variants_quantitative_scores: {quantitative_score: {name: {_eq: "SIFT"}, value: {_lt: 0.05}}}
      neighbor_edits: {is_best: {_eq: true}, improves: {_eq: true}}
"""

SETS = {
    "treatable": {"where": WHERE_TREATABLE, "prefix": "TableS1_Fixable_variants",
                  "xlsx": "TableS1_Fix_Supplementary_Tables.xlsx", "extra": None},
    "rescue": {"where": WHERE_RESCUE, "prefix": "TableS2_Rescue_variants",
               "xlsx": "TableS2_Rescue_Supplementary_Tables.xlsx", "extra": "rescue"},
    "improve": {"where": WHERE_IMPROVE, "prefix": "TableS3_Missense_optimization_variants",
                "xlsx": "TableS3_Missense_Optimization_Supplementary_Tables.xlsx", "extra": "improve"},
}

QUERY_TEMPLATE = """
query VariantTable($limit: Int!, $offset: Int!) {
  result: variants(
    where: {
__WHERE__
    }
    order_by: {id: asc}
    limit: $limit
    offset: $offset
  ) {
    id
    coordinate { chr start end }
    ref
    alt
    ref_codon { nt1 nt2 nt3 }
    alt_codon { nt1 nt2 nt3 }
    variants_consequences { consequence { name } }
    gene {
      ensg
      symbol
      strand
      genes_quantitative_scores(where: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}) {
        quantitative_score { value }
      }
      genes_qualitative_scores(where: {qualitative_score: {name: {_eq: "SFARI Syndromic"}}}) {
        qualitative_score { value }
      }
    }
    cds_position
    variants_quantitative_scores(where: {quantitative_score: {name: {_in: ["SIFT", "CADD_phred"]}}}) {
      quantitative_score { name value }
    }
    allele_frequency { gnomade_af max_af max_af_pops }
    best_edit: neighbor_edits(where: {is_best: {_eq: true}}) {
      edited_positions
      restores_reference
      pre_codon { nt1 nt2 nt3 amino_acid { short_name } }
      post_codon { nt1 nt2 nt3 amino_acid { short_name } }
      neighbor_edit_scores(where: {source: {_eq: "SIFT"}}) { value improvement }
    }
    variants_guides(where: {guide: {type: {_eq: "pre_mRNA"}}}) {
      guide { hits_85 }
      rescue_sift
      bystanders {
        id
        genomic_position
        bystanders_consequences { consequence { name } }
        cds_position
        in_cds
        sift_score
        sift_prediction
        cadd_phred
        ref_codon { nt1 nt2 nt3 amino_acid { short_name } }
        alt_codon { nt1 nt2 nt3 amino_acid { short_name } }
      }
    }
  }
}
"""

# Consequence severity ranking, most-severe first. StopGained outranks Splice, matching the
# rule the analysis actually applies (`New_Consequence` in
# figures/ver1/00_prep_data.rmd, which tests StopGained first and is inherited
# by every figure). See docs/decisions/0001-consequence-precedence.md.
#
# This previously listed Splice first, mirroring a `variants_rank` constant in pipeline/phase2_db.py
# that was itself never read and has been removed. No published number depended on either:
# applied here it collapses *bystander* consequences, and no bystander in the dataset carries
# both Splice and StopGained.
CONSEQUENCE_RANK = [
    "StopGained", "Splice", "Frameshift", "StopLost", "StartLost", "Missense",
    "Synonymous", "3UTR", "5UTR", "Intron", "Upstream", "Downstream", "Other",
]
_RANK_INDEX = {name: i for i, name in enumerate(CONSEQUENCE_RANK)}

# Bystander tuple (merged + wide): data keys vs. display names (position base noted).
BYSTANDER_TUPLE_FIELDS = [
    "consequence", "genomic_position", "sift_score", "sift_prediction",
    "cadd_phred", "ref_codon", "ref_aa", "alt_codon", "alt_aa",
]
BYSTANDER_TUPLE_HEADER = [
    "consequence", "genomic_position (1-based)", "sift_score", "sift_prediction",
    "cadd_phred", "ref_codon", "ref_aa", "alt_codon", "alt_aa",
]

# A +-10 nt editing window has at most 20 possible bystander A positions, so the
# wide CSV always uses 20 fixed bystander columns (unused slots are NULL).
WINDOW_MAX_BYSTANDERS = 20

# Per-variant scalar columns, shared in the SAME order by all three formats.
VARIANT_BASE_COLUMNS = [
    "variant_id", "chr", "pos (1-based)", "ref", "alt", "ref_codon", "alt_codon",
    "consequences", "gene_ensg", "gene_symbol", "strand", "sfari_gene_score",
    "sfari_syndromic", "sift", "cadd", "gnomade_af", "max_af", "max_af_pops", "off-targets",
]

# Extra columns only on the rescue set.
RESCUE_EXTRA_COLUMNS = ["rescue_sift"]

# Extra columns only on the missense-optimization set — describe the chosen (best) neighbor edit:
# the variant's own SIFT is already in the shared `sift` column; these add the edit outcome.
# `alt_aa` is the patient (alt) codon's amino acid.
IMPROVE_EXTRA_COLUMNS = [
    "edited_codon_positions", "alt_aa", "edited_codon",
    "edited_aa", "edit_sift", "sift_gain",
]

# Map a set's `extra` key to its extra column list (None -> no extra columns).
EXTRA_COLUMNS = {"rescue": RESCUE_EXTRA_COLUMNS, "improve": IMPROVE_EXTRA_COLUMNS}

def variant_columns(extra):
    return VARIANT_BASE_COLUMNS + EXTRA_COLUMNS.get(extra, [])


def most_severe(names):
    """Pick the most-severe consequence by CONSEQUENCE_RANK; unknown terms sort last."""
    names = [n for n in names if n]
    if not names:
        return None
    return min(names, key=lambda n: _RANK_INDEX.get(n, len(CONSEQUENCE_RANK)))


def _codon_str(codon):
    """Concatenate nt1+nt2+nt3 into a codon string, or None if absent."""
    if not codon:
        return None
    nts = [codon.get("nt1"), codon.get("nt2"), codon.get("nt3")]
    if all(n is None for n in nts):
        return None
    return "".join(n or "" for n in nts)


def _codon_aa(codon):
    if not codon:
        return None
    aa = codon.get("amino_acid")
    return aa.get("short_name") if aa else None


def _first_score(rows, container_key):
    """First .{container_key}.value from a list of score join rows (or None)."""
    for row in rows or []:
        sc = row.get(container_key)
        if sc and sc.get("value") is not None:
            return sc["value"]
    return None


def _first(values):
    return values[0] if values else None


NULL = "NULL"  # explicit placeholder for empty values (no blank cells anywhere)


def _fmt(x):
    """Render a scalar for a CSV cell: empty -> NULL, bool -> True/False, else str."""
    if isinstance(x, bool):
        return "True" if x else "False"
    if x is None or x == "":
        return NULL
    return str(x)


def _join_scalars(values):
    """';'-join scalars (comma-free so the cell is never quoted); empty list -> NULL."""
    if not values:
        return NULL
    return ";".join(_fmt(x) for x in values)


def _num(x):
    """Coerce a numeric-looking value to int/float for xlsx (so Excel sorts it as
    a number); leave genuinely non-numeric values (e.g. chr 'X') as-is."""
    if x is None or isinstance(x, (int, float)):
        return x
    s = str(x).strip()
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return x


def _nullify(x):
    """xlsx cell: empty -> 'NULL'; booleans and numbers are kept as-is."""
    if isinstance(x, bool):
        return x
    if x is None or x == "":
        return NULL
    return x


def extract(v):
    """Normalize one GraphQL variant into a flat dict + list of bystander dicts."""
    gene = v.get("gene") or {}
    coord = v.get("coordinate") or {}

    score_by_name = {}
    for row in v.get("variants_quantitative_scores") or []:
        sc = row.get("quantitative_score") or {}
        if sc.get("name") is not None:
            score_by_name[sc["name"]] = sc.get("value")

    consequences = [
        (c.get("consequence") or {}).get("name")
        for c in (v.get("variants_consequences") or [])
    ]

    # pre_mRNA guide only (query-filtered): off-targets, rescue_sift, bystanders.
    guides = v.get("variants_guides") or []
    off_targets = [(g.get("guide") or {}).get("hits_85") for g in guides]
    rescue_sift = [g.get("rescue_sift") for g in guides]

    bystanders = []
    for g in guides:
        for b in g.get("bystanders") or []:
            bystanders.append({
                "bystander_id": b.get("id"),
                "genomic_position": b.get("genomic_position"),
                "consequence": most_severe(
                    (bc.get("consequence") or {}).get("name")
                    for bc in (b.get("bystanders_consequences") or [])
                ),
                "cds_position": b.get("cds_position"),
                "in_cds": b.get("in_cds"),
                "sift_score": b.get("sift_score"),
                "sift_prediction": b.get("sift_prediction"),
                "cadd_phred": b.get("cadd_phred"),
                "ref_codon": _codon_str(b.get("ref_codon")),
                "ref_aa": _codon_aa(b.get("ref_codon")),
                "alt_codon": _codon_str(b.get("alt_codon")),
                "alt_aa": _codon_aa(b.get("alt_codon")),
            })

    # best neighbor-drill edit (Improve set only); empty for treatable/rescue variants.
    be = _first(v.get("best_edit")) or {}
    be_nes = _first(be.get("neighbor_edit_scores")) or {}

    return {
        "variant_id": v.get("id"),
        "chr": coord.get("chr"),
        "pos": coord.get("start"),  # start == end for SNVs; 1-based (verified vs hg19)
        "ref": v.get("ref"),
        "alt": v.get("alt"),
        "ref_codon": _codon_str(v.get("ref_codon")),
        "alt_codon": _codon_str(v.get("alt_codon")),
        "consequences": consequences,
        "gene_ensg": gene.get("ensg"),
        "gene_symbol": gene.get("symbol"),
        "strand": gene.get("strand"),
        "sfari_gene_score": _first_score(gene.get("genes_quantitative_scores"), "quantitative_score"),
        "sfari_syndromic": _first_score(gene.get("genes_qualitative_scores"), "qualitative_score"),
        "cds_position": v.get("cds_position"),
        "sift": score_by_name.get("SIFT"),
        "cadd": score_by_name.get("CADD_phred"),
        "gnomade_af": (v.get("allele_frequency") or {}).get("gnomade_af"),
        "max_af": (v.get("allele_frequency") or {}).get("max_af"),
        "max_af_pops": (v.get("allele_frequency") or {}).get("max_af_pops"),
        "off_targets": off_targets,
        "rescue_sift": rescue_sift,
        "improve_edited_codon_positions": be.get("edited_positions"),
        "improve_patient_aa": _codon_aa(be.get("pre_codon")),
        "improve_edited_codon": _codon_str(be.get("post_codon")),
        "improve_edited_aa": _codon_aa(be.get("post_codon")),
        "improve_edit_sift": be_nes.get("value"),
        "improve_sift_gain": be_nes.get("improvement"),
        "bystanders": bystanders,
    }


def _bystander_tuple(b):
    """'(consequence;cds_position;...)' — comma-free, never quoted by csv."""
    return "(" + ";".join(_fmt(b[k]) for k in BYSTANDER_TUPLE_FIELDS) + ")"


def _base_row(d, extra):
    """Per-variant scalar cells (same order as VARIANT_BASE_COLUMNS) for CSV output."""
    row = {
        "variant_id": _fmt(d["variant_id"]), "chr": _fmt(d["chr"]), "pos (1-based)": _fmt(d["pos"]),
        "ref": _fmt(d["ref"]), "alt": _fmt(d["alt"]),
        "ref_codon": _fmt(d["ref_codon"]), "alt_codon": _fmt(d["alt_codon"]),
        "consequences": _join_scalars(d["consequences"]),
        "gene_ensg": _fmt(d["gene_ensg"]), "gene_symbol": _fmt(d["gene_symbol"]), "strand": _fmt(d["strand"]),
        "sfari_gene_score": _fmt(d["sfari_gene_score"]), "sfari_syndromic": _fmt(d["sfari_syndromic"]),
        "sift": _fmt(d["sift"]), "cadd": _fmt(d["cadd"]),
        "gnomade_af": _fmt(d["gnomade_af"]),
        "max_af": _fmt(d["max_af"]), "max_af_pops": _fmt(d["max_af_pops"]),
        "off-targets": _join_scalars(d["off_targets"]),
    }
    if extra == "rescue":
        row["rescue_sift"] = _join_scalars(d["rescue_sift"])
    elif extra == "improve":
        row["edited_codon_positions"] = _join_scalars(d["improve_edited_codon_positions"])
        row["alt_aa"] = _fmt(d["improve_patient_aa"])
        row["edited_codon"] = _fmt(d["improve_edited_codon"])
        row["edited_aa"] = _fmt(d["improve_edited_aa"])
        row["edit_sift"] = _fmt(d["improve_edit_sift"])
        row["sift_gain"] = _fmt(d["improve_sift_gain"])
    return row


def write_merged_csv(path, data, extra):
    by_col = f"bystanders({'; '.join(BYSTANDER_TUPLE_HEADER)})"
    cols = variant_columns(extra) + ["bystanders_count", by_col]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for d in data:
            row = _base_row(d, extra)
            row["bystanders_count"] = len(d["bystanders"])
            row[by_col] = " | ".join(_bystander_tuple(b) for b in d["bystanders"]) or NULL
            w.writerow(row)
    logging.info(f"  ✅ MERGED  -> {path}  ({len(data)} variant rows)")


def write_wide_csv(path, data, max_bystanders, extra):
    # Fixed 20 bystander columns (the +-10 nt window's theoretical max); unused are
    # NULL. Extend only if data exceeds 20, so nothing is ever silently dropped.
    n_cols = max(WINDOW_MAX_BYSTANDERS, max_bystanders)
    if max_bystanders > WINDOW_MAX_BYSTANDERS:
        logging.warning(f"  a variant has {max_bystanders} bystanders, exceeding the "
                        f"+-10 nt max of {WINDOW_MAX_BYSTANDERS}; using {n_cols} columns.")
    first_hdr = f"bystander_1 ({'; '.join(BYSTANDER_TUPLE_HEADER)})"
    by_cols = [first_hdr] + [f"bystander_{i}" for i in range(2, n_cols + 1)]
    cols = variant_columns(extra) + ["bystanders_count"] + by_cols
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for d in data:
            row = _base_row(d, extra)
            row["bystanders_count"] = len(d["bystanders"])
            tuples = [_bystander_tuple(b) for b in d["bystanders"]]
            for i, col in enumerate(by_cols):
                row[col] = tuples[i] if i < len(tuples) else "NULL"
            w.writerow(row)
    logging.info(f"  ✅ WIDE    -> {path}  ({len(data)} rows, {len(by_cols)} bystander cols, NULL-filled)")


def write_xlsx(path, data, extra):
    """Single-sheet workbook mirroring the merged CSV: one row per variant, all
    bystanders packed into one column. xlsx carries no delimiters, so it renders
    cleanly in spreadsheet apps (Google Sheets/Drive, Excel) — no quoting needed."""
    by_col = f"bystanders({'; '.join(BYSTANDER_TUPLE_HEADER)})"
    cols = variant_columns(extra) + ["bystanders_count", by_col]
    wb = Workbook()
    ws = wb.active
    ws.title = "Variants"
    ws.append(cols)
    for d in data:
        row = _base_row(d, extra)
        row["bystanders_count"] = len(d["bystanders"])
        row[by_col] = " | ".join(_bystander_tuple(b) for b in d["bystanders"]) or NULL
        ws.append([_num(row[c]) for c in cols])

    # Light publication styling: bold header, frozen header row, autofilter, widths.
    for cell in ws[1]:
        cell.font = Font(bold=True)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{ws.max_row}"
    for i, name in enumerate(cols, 1):
        ws.column_dimensions[get_column_letter(i)].width = min(max(len(name) + 2, 11), 40)

    wb.save(path)
    logging.info(f"  ✅ XLSX    -> {path}  ({len(data)} rows, single sheet)")


def fetch_variants(query, limit_total, page_size):
    rows, offset = [], 0
    while True:
        page = page_size if limit_total is None else min(page_size, limit_total - len(rows))
        if page <= 0:
            break
        batch = gql(query, {"limit": page, "offset": offset})["result"]
        rows.extend(batch)
        logging.info(f"  fetched {len(rows)} variants...")
        if len(batch) < page:
            break
        offset += page
    return rows


def run_set(set_name, fmt, out_dir, limit, page_size):
    cfg = SETS[set_name]
    logging.info(f"=== set '{set_name}' ===")
    query = QUERY_TEMPLATE.replace("__WHERE__", cfg["where"])
    variants = fetch_variants(query, limit, page_size)
    data = [extract(v) for v in variants]
    max_bystanders = max((len(d["bystanders"]) for d in data), default=0)
    logging.info(f"  {len(data)} variants; max bystanders/variant: {max_bystanders}")

    extra = cfg["extra"]
    if fmt in ("all", "merged"):
        write_merged_csv(out_dir / f"{cfg['prefix']}.csv", data, extra)
    if fmt in ("all", "wide"):
        write_wide_csv(out_dir / f"{cfg['prefix']}_wide.csv", data, max_bystanders, extra)
    if fmt in ("all", "xlsx"):
        write_xlsx(out_dir / cfg["xlsx"], data, extra)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=["all", "treatable", "rescue", "improve"], default="all",
                    help="Which variant set(s) to export (default: all)")
    ap.add_argument("--format", choices=["all", "merged", "wide", "xlsx"], default="merged",
                    help="Which output format(s) to write (default: merged CSV only)")
    ap.add_argument("--out-dir", default=str(PROJECT_ROOT / "Output"),
                    help="Output directory (default: Output/)")
    ap.add_argument("--limit", type=int, default=None, help="Max variants per set (default: all)")
    ap.add_argument("--page-size", type=int, default=200, help="Variants per GraphQL request")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sets = ["treatable", "rescue", "improve"] if args.set == "all" else [args.set]
    for set_name in sets:
        run_set(set_name, args.format, out_dir, args.limit, args.page_size)


if __name__ == "__main__":
    main()
