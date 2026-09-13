"""
Populate the `neighbor_edits` table for non-G>A missense variants, scored fully offline (SIFT,
release 113 — same version as the variants' stored SIFT).

Per variant: enumerate A→G codon edits (core) → map each codon's CDS positions to genomic via the
release-87 GTF and express each edit as a genomic delins (verified against the variant's own
coordinate + the genome FASTA) → offline VEP for SIFT → score, pick best → insert.

Variants we can't verify (mapper mismatch / split codon / FASTA mismatch) or that have no SIFT are
recorded WITHOUT a score. Idempotent: variants that already have neighbor_edits are excluded.

`populate_neighbor_edits()` is called by phase3 (reproducible) and by populate_now (backfill); it
reads paths from the already-loaded environment.
"""

import os
import sys
import json
import logging
from pathlib import Path

from pyfaidx import Fasta
from tqdm import tqdm

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
import helpers as hp
from . import core as nd
from . import offline_vep as ov

log = logging.getLogger(__name__)


# _ROOT is pipeline/ - right for the sys.path insert above, since helpers.py lives there.
# Data paths are relative to the REPOSITORY, one level further up. Conflating the two
# resolved Resources/... to pipeline/Resources/... once neighbor_drill moved a level deeper.
_REPO = Path(__file__).resolve().parents[2]


def _resolve(p):
    p = Path(p).expanduser()
    return p if p.is_absolute() else (_REPO / p).resolve()


GET_IMPROVE_VARIANTS = """
query GetImproveVariants($limit: Int!, $offset: Int!) {
  variants(
    where: {
      class: {_eq: "SNV"}
      variants_consequences: {consequence: {name: {_eq: "Missense"}}}
      alt_codon: {_or: [{nt1: {_eq: "A"}}, {nt2: {_eq: "A"}}, {nt3: {_eq: "A"}}]}
      _not: {_or: [
        {_and: [{ref: {_eq: "G"}}, {alt: {_eq: "A"}}, {_or: [{gene: {strand: {_eq: "+"}}}, {gene: {strand: {_is_null: true}}}]}]}
        {_and: [{ref: {_eq: "C"}}, {alt: {_eq: "T"}}, {gene: {strand: {_eq: "-"}}}]}
        {neighbor_edits: {}}
      ]}
    }
    order_by: {id: asc}
    limit: $limit
    offset: $offset
  ) {
    id cds_position position_in_codon
    coordinate { chr start }
    gene { strand }
    variants_features(where: {feature: {type: {_eq: "Transcript"}}}) { feature { identifier version_number } }
    ref_codon { nt1 nt2 nt3 amino_acid { short_name } }
    alt_codon { id nt1 nt2 nt3 }
    sift: variants_quantitative_scores(where: {quantitative_score: {name: {_eq: "SIFT"}}}) { quantitative_score { value } }
  }
}
"""

INSERT = """
mutation InsertNeighborEdits($objects: [neighbor_edits_insert_input!]!) {
  insert_neighbor_edits(objects: $objects,
    on_conflict: {constraint: neighbor_edits_variant_id_edited_positions_key, update_columns: []}
  ) { affected_rows }
}
"""


def _fetch_variants(gql, limit):
    rows, offset, page = [], 0, 1000
    while True:
        ps = page if limit is None else min(page, limit - len(rows))
        if ps <= 0:
            break
        batch = gql(GET_IMPROVE_VARIANTS, {"limit": ps, "offset": offset})["variants"]
        rows.extend(batch)
        if len(batch) < ps:
            break
        offset += ps
    return rows


def populate_neighbor_edits(limit=None, dry_run=False, gtf=None, forks=8, insert_batch=500):
    """Compute + insert neighbor_edits. Reads GENOME_REFERENCE_HG19 / VEP_CACHE / OUTPUT_DIR from env."""
    gql = hp.gql
    scorer = nd.SiftScorer()
    genome = _resolve(os.environ["GENOME_REFERENCE_HG19"])
    vep_cache = _resolve(os.environ["VEP_CACHE"])
    work = _resolve(os.environ["OUTPUT_DIR"]) / "neighbor_drill"
    work.mkdir(parents=True, exist_ok=True)
    gtf = _resolve(gtf or os.environ.get("NEIGHBOR_DRILL_GTF", "Resources/Homo_sapiens.GRCh37.87.gtf.gz"))

    codon_lookup = {ov.to_dna(f"{c['nt1']}{c['nt2']}{c['nt3']}"): c["id"]
                    for c in gql("query { codons { id nt1 nt2 nt3 } }", {})["codons"]}
    variants = _fetch_variants(gql, limit)
    log.info(f"[neighbor_edits] candidate variants: {len(variants)}")
    if not variants:
        return

    fa = Fasta(str(genome), as_raw=True)
    chr_prefix = "chr" if any(k.startswith("chr") for k in fa.keys()) else ""
    wanted = {v["variants_features"][0]["feature"]["identifier"]
              for v in variants if v.get("variants_features")}
    tx = ov.load_cds_intervals(gtf, wanted)
    log.info(f"[neighbor_edits] loaded CDS for {len(tx)}/{len(wanted)} transcripts")

    # ── Phase 1: build per-option genomic delins ──
    plans, vcf_rows = [], []
    skip = {"no_tx": 0, "guard": 0, "split": 0, "fasta": 0, "no_options": 0, "bad_codon": 0}
    for v in variants:
        feats = v.get("variants_features") or []
        feat = (feats[0].get("feature") if feats else None) or {}
        enst = feat.get("identifier")
        t = tx.get(enst)
        if t is None or v.get("cds_position") is None or v.get("position_in_codon") is None:
            skip["no_tx"] += 1
            continue
        cds_pos, pic = int(v["cds_position"]), int(v["position_in_codon"])
        chrom = str(v["coordinate"]["chr"])
        if ov.cds_to_genomic(t, cds_pos) != v["coordinate"]["start"]:
            skip["guard"] += 1
            continue
        cds_start = cds_pos - (pic - 1)
        gs = [ov.cds_to_genomic(t, cds_start + i) for i in range(3)]
        if any(g is None for g in gs):
            skip["bad_codon"] += 1
            continue
        strand = t["strand"]
        if strand == "+":
            lo, contiguous = gs[0], (gs[1] == gs[0] + 1 and gs[2] == gs[0] + 2)
        else:
            lo, contiguous = gs[2], (gs[1] == gs[0] - 1 and gs[2] == gs[0] - 2)
        if not contiguous:
            skip["split"] += 1
            continue

        ref_codon = ov.to_dna(f"{v['ref_codon']['nt1']}{v['ref_codon']['nt2']}{v['ref_codon']['nt3']}")
        alt_codon = ov.to_dna(f"{v['alt_codon']['nt1']}{v['alt_codon']['nt2']}{v['alt_codon']['nt3']}")
        ref_aa = (v["ref_codon"].get("amino_acid") or {}).get("short_name")
        genomic_ref = ref_codon if strand == "+" else ov.rc(ref_codon)
        if str(fa[chr_prefix + chrom][lo - 1: lo + 2]).upper() != genomic_ref:
            skip["fasta"] += 1
            continue

        options = nd.enumerate_neighbor_edits(alt_codon, ref_aa)
        if not options:
            skip["no_options"] += 1
            continue

        baseline = v["sift"][0]["quantitative_score"]["value"] if v.get("sift") else None
        plan = {"vid": v["id"], "enst": enst, "pic": pic, "strand": strand,
                "pre_codon_id": v["alt_codon"]["id"], "baseline": baseline,
                "options": options, "post_ids": [], "edited_cds": []}
        for i, o in enumerate(options):
            post_dna = ov.to_dna(o.post_codon)
            genomic_alt = post_dna if strand == "+" else ov.rc(post_dna)
            vcf_rows.append((chrom, lo, f"ne_{v['id']}_{i}", genomic_ref, genomic_alt))
            plan["post_ids"].append(codon_lookup.get(post_dna))
            plan["edited_cds"].append([cds_start + (p - 1) for p in o.edited_positions])
        plans.append(plan)

    log.info(f"[neighbor_edits] plannable: {len(plans)} | options: {len(vcf_rows)} | skipped: {skip}")
    if not plans:
        return

    # ── Phase 2: write VCF + offline VEP ──
    vcf_path, out_path = work / "neighbor_edits.vcf", work / "neighbor_edits.vep.txt"
    vcf_rows.sort(key=lambda r: (r[0], r[1]))
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, vid, ref, alt in vcf_rows:
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")
    log.info(f"[neighbor_edits] running offline VEP on {len(vcf_rows)} records ...")
    ov.run_offline_vep(vcf_path, out_path, vep_cache=vep_cache, genome_fasta=genome, work_dir=work, forks=forks)
    vep = ov.parse_vep(out_path)

    # ── Phase 3: score + build insert objects ──
    objects, stat = [], {"scored": 0, "restorers": 0, "improves": 0, "is_best": 0, "no_sift_var": 0}
    for plan in plans:
        regions_per, any_sift = [], False
        for i, o in enumerate(plan["options"]):
            row = ov.row_for_transcript(vep.get(f"ne_{plan['vid']}_{i}", []), plan["enst"])
            label, score = ov.parse_sift(row.get("Extra", "")) if row else (None, None)
            nd.set_score(o, scorer, score, label, plan["baseline"])
            regions_per.append(hp.consequence_terms_to_regions((row.get("Consequence", "") or "").split(",") if row else []))
            any_sift = any_sift or (score is not None)
        nd.evaluate_options(plan["options"], scorer, plan["pic"])
        if not any_sift and not any(o.restores_reference for o in plan["options"]):
            stat["no_sift_var"] += 1
        for i, o in enumerate(plan["options"]):
            stat["scored"] += 1 if (o.scores.get("SIFT") and o.scores["SIFT"].value is not None) else 0
            stat["restorers"] += o.restores_reference
            stat["improves"] += o.improves
            stat["is_best"] += o.is_best
            obj = {"variant_id": plan["vid"], "pre_codon_id": plan["pre_codon_id"],
                   "post_codon_id": plan["post_ids"][i], "edited_positions": list(o.edited_positions),
                   "edited_cds_positions": plan["edited_cds"][i], "restores_reference": o.restores_reference,
                   "is_best": o.is_best, "improves": o.improves}
            sc = [{"source": s.source, "value": s.value, "label": s.label, "improvement": s.improvement}
                  for s in o.scores.values() if s.value is not None]
            if sc:
                obj["neighbor_edit_scores"] = {"data": sc}
            if regions_per[i]:
                obj["neighbor_edits_consequences"] = {"data": [
                    {"consequence": {"data": {"name": r},
                                     "on_conflict": {"constraint": "consequences_name_key", "update_columns": ["name"]}}}
                    for r in regions_per[i]]}
            objects.append(obj)

    log.info(f"[neighbor_edits] built {len(objects)} option rows | {stat}")
    if dry_run:
        log.info("[neighbor_edits] --dry-run: not writing.")
        return

    total = 0
    for i in tqdm(range(0, len(objects), insert_batch), desc="Insert neighbor_edits"):
        total += gql(INSERT, {"objects": objects[i:i + insert_batch]})["insert_neighbor_edits"]["affected_rows"]
    log.info(f"[neighbor_edits] ✅ inserted {total} rows.")
