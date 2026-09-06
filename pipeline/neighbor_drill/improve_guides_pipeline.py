"""
Build guides + bystanders for the "Improve" editing class — SFARI protein-coding non-G>A missense
variants that have a scored best neighbor-drill edit (is_best).

Scope + scoring (locked spec):
  * Variant set: SFARI protein-coding, non-G>A, Missense, with an is_best neighbor edit (~1,923).
  * Guides + BLAT off-target: built around the VARIANT position (same as Fix/Rescue).
  * Bystanders: ±10 nt editable A's, EXCLUDING the best edit's intended edit positions (those are
    the therapy, not bystander events).
  * Bystander scoring is ONLINE (Ensembl GRCh37 VEP REST, ?CADD=1&pick=1), exactly like Fix/Rescue:
      - CADD (cadd_raw/cadd_phred): the lone A→G SNV — nucleotide level.
      - SIFT (sift_score/sift_prediction): the amino-acid consequence of the lone A→G — EXCEPT for
        the ≤1 bystander that sits in the variant's own codon, whose SIFT is computed codon-aware
        via a coding delins on the therapy codon (post_codon) + the bystander edit (e.g. GAA→CGG).
      - AlphaMissense is not requested/stored (am_* left null).
  * Marginal/net SIFT is left to the analyst — derivable from the already-stored variant SIFT and
    the intended-edit SIFT (neighbor_edit_scores); no extra columns.

The guide-extraction / bystander-finding / BLAT helpers are copied from phase3 (self-contained).
`populate_improve_guides()` is called by phase3 (reproducible) and by populate_now (backfill).
"""

import os
import sys
import json
import time
import math
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import pandas as pd
from Bio.Seq import Seq
from pyfaidx import Fasta
from tqdm import tqdm

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
import helpers as hp
from . import offline_vep as ov

log = logging.getLogger(__name__)

# module-level config, set by populate_improve_guides()
genome_fasta = cds_fasta = None
needs_chr_prefix = False
BLAT_PATH = None
BLAT_WORKERS = 8
BLAT_DIR = None
VEP_SERVER = None
VEP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
VEP_HGVS_EXT = "/vep/human/hgvs?CADD=1&pick=1"   # CADD only; no AlphaMissense (per spec)

BLAT_CHUNK_SIZE = 500
IDENTITY_THRESHOLDS = [100, 95, 90, 85]
col_names = ['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert',
             'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd',
             'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts']


def _resolve(p):
    p = Path(p).expanduser()
    return p if p.is_absolute() else (_ROOT / p).resolve()


# ── guide extraction + bystander finding (verbatim from phase3) ──
def extract_guide(fasta_handle, key, pos, strand="+", n=20):
    full_seq = fasta_handle[key]
    seq_len = len(full_seq)
    idx = pos - 1
    if idx < 0 or idx >= seq_len:
        raise ValueError(f"Position {pos} out of bounds for {key} ({seq_len})")
    start = max(0, idx - n)
    end = min(seq_len, idx + n + 1)
    guide_seq = str(full_seq[start:end])
    variant_idx = idx - start
    if strand == "-":
        guide_seq = str(Seq(guide_seq).reverse_complement())
        variant_idx = len(guide_seq) - 1 - variant_idx
    return guide_seq.upper(), variant_idx


def find_bystanders(guide, variant_idx, position, strand="+", window_size=10):
    """Editable A's within ±window_size of the variant (variant itself excluded)."""
    position = int(position)
    context = guide.upper().replace("T", "U")
    context_len = len(context)
    out = []
    for i in range(max(0, variant_idx - window_size), min(context_len, variant_idx + window_size + 1)):
        if i == variant_idx:
            continue
        if context[i] == "A":
            offset = i - variant_idx
            if strand == "-":
                offset = -offset
            out.append({"has_u_before": (i > 0 and context[i - 1] == "U"),
                        "has_u_after": (i < context_len - 1 and context[i + 1] == "U"),
                        "position": position + offset})
    return out


def make_guide_obj(guide_type, variant_id, sequence, variant_idx, edit_type=None):
    return {"hits_85": None, "hits_90": None, "hits_95": None, "hits_100": None,
            "sequence": sequence, "variant_id": variant_id, "variant_idx": variant_idx,
            "bystanders": [], "type": guide_type, "edit_type": edit_type, "rescue_sift": None}


def process_variant_cpu(v, edit_type="Improve"):
    """Guides + bystanders around the VARIANT (no network)."""
    vid = v["id"]
    strand = (v["gene"] or {}).get("strand") or "+"
    chrom = str(v["coordinate"]["chr"])
    if needs_chr_prefix and not chrom.startswith("chr"):
        chrom = "chr" + chrom
    gpos = v["coordinate"]["start"]
    res = {"variant_id": vid, "strand": strand, "chrom": chrom, "transcript_id": None, "tid_ver": None,
           "pre_bystanders": [], "rna_bystanders": [], "dna_guide_obj": None, "rna_guide_obj": None}
    try:
        dna_guide, dna_vi = extract_guide(genome_fasta, chrom, gpos, strand, 20)
    except Exception as e:
        res["error"] = f"pre_mRNA: {e}"
        return res
    res["dna_guide_obj"] = make_guide_obj("pre_mRNA", vid, dna_guide, dna_vi + 1, edit_type)
    res["pre_bystanders"] = find_bystanders(dna_guide, dna_vi, gpos, strand)
    feats = v.get("variants_features") or []
    feat = (feats[0].get("feature") if feats else None) or {}
    tid = feat.get("identifier")
    if tid and v.get("cds_position"):
        tid_ver = f"{tid}.{feat['version_number']}" if feat.get("version_number") else tid
        try:
            cds_pos = int(v["cds_position"])
            rna_guide, rna_vi = extract_guide(cds_fasta, tid_ver, cds_pos)
            res["rna_guide_obj"] = make_guide_obj("mature_mRNA", vid, rna_guide, rna_vi + 1, edit_type)
            res["rna_bystanders"] = find_bystanders(rna_guide, rna_vi, cds_pos)
            res["transcript_id"] = tid
            res["tid_ver"] = tid_ver
        except Exception as e:
            res["error_rna"] = f"mature: {e}"
    return res


# ── BLAT (from phase3; ThreadPool since BLAT is a subprocess) ──
def run_blat_chunk(chunk_guides, reference, chunk_id, guide_type):
    query_path = BLAT_DIR / f"improve_{guide_type}_chunk{chunk_id}.fa"
    results_path = BLAT_DIR / f"improve_{guide_type}_chunk{chunk_id}.psl"
    with open(query_path, "w") as f:
        for g in chunk_guides:
            f.write(f">{g['variant_id']}\n{g['sequence']}\n")
    import subprocess
    subprocess.run([BLAT_PATH, str(reference), str(query_path), str(results_path),
                    "-out=psl", "-stepSize=5", "-repMatch=2253", "-minScore=20", "-minIdentity=0"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        df = pd.read_csv(results_path, sep="\t", skiprows=5, names=col_names)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=col_names)
    os.remove(query_path)
    os.remove(results_path)
    return df


def blat_all_guides(guides, reference, guide_type):
    if not guides:
        return {}
    chunks = [guides[i:i + BLAT_CHUNK_SIZE] for i in range(0, len(guides), BLAT_CHUNK_SIZE)]
    all_dfs = []
    with ThreadPoolExecutor(max_workers=min(BLAT_WORKERS, len(chunks))) as pool:
        futures = {pool.submit(run_blat_chunk, ch, reference, idx, guide_type): idx for idx, ch in enumerate(chunks)}
        for fut in tqdm(as_completed(futures), total=len(futures), desc=f"BLAT {guide_type}"):
            all_dfs.append(fut.result())
    blat_df = pd.concat(all_dfs, ignore_index=True)
    if blat_df.empty:
        return {}
    blat_df["identity"] = blat_df.apply(lambda r: hp.calculate_psl_percent_identity(
        r["matches"], r["misMatches"], r["repMatches"], r["qNumInsert"], r["tNumInsert"],
        r["qStart"], r["qEnd"], r["tStart"], r["tEnd"]), axis=1)
    blat_df["qName"] = blat_df["qName"].astype(str)
    hits = {}
    for thr in IDENTITY_THRESHOLDS:
        for vid, n in blat_df[blat_df["identity"] >= thr]["qName"].value_counts().to_dict().items():
            hits.setdefault(vid, {})[f"hits_{thr}"] = max(0, n - 1)
    for vid in hits:
        for thr in IDENTITY_THRESHOLDS:
            hits[vid].setdefault(f"hits_{thr}", 0)
    return hits


# ── online VEP (HGVS batch, CADD+SIFT) ──
def _vep_hgvs(hgvs_list, batch_size=200, max_workers=8):
    """POST HGVS in batches to GRCh37 VEP; return {input_hgvs: entry}."""
    uniq = sorted(set(hgvs_list))
    out = {}
    starts = list(range(0, len(uniq), batch_size))

    def _send(start):
        chunk = uniq[start:start + batch_size]
        payload = json.dumps({"hgvs_notations": chunk})
        for attempt in range(5):
            r = requests.post(VEP_SERVER + VEP_HGVS_EXT, headers=VEP_HEADERS, data=payload, timeout=120)
            if r.status_code == 200:
                return {e.get("input"): e for e in r.json()}
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 1))); continue
            if r.status_code in (502, 503, 504):
                time.sleep(0.5 * (2 ** attempt)); continue
            raise RuntimeError(f"VEP {r.status_code}: {r.text[:500]}")
        raise RuntimeError(f"VEP failed after retries (batch at {start})")

    with ThreadPoolExecutor(max_workers=min(max_workers, max(1, len(starts)))) as pool:
        futures = [pool.submit(_send, s) for s in starts]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="VEP bystanders"):
            out.update(fut.result())
    return out


def _pick_tc(entry, transcript_id=None):
    if entry is None:
        return {}
    tcs = entry.get("transcript_consequences", []) or []
    if transcript_id:
        m = next((t for t in tcs if t.get("transcript_id") == transcript_id), None)
        if m:
            return m
    return tcs[0] if tcs else {}


GET_IMPROVE_GUIDE_VARIANTS = """
query($limit: Int!, $offset: Int!) {
  variants(
    where: {
      class: {_eq: "SNV"}
      _not: {_or: [
        {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
        {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
        {variants_guides: {edit_type: {_eq: "Improve"}}}
      ]}
      variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}}
      gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}}
      variants_consequences: {consequence: {name: {_eq: "Missense"}}}
      neighbor_edits: {is_best: {_eq: true}}
    }
    order_by: {id: asc} limit: $limit offset: $offset
  ) {
    id cds_position position_in_codon
    gene { strand }
    coordinate { chr start }
    variants_features(where: {feature: {type: {_eq: "Transcript"}}}) { feature { identifier version_number } }
    ref_codon { nt1 nt2 nt3 }
    best: neighbor_edits(where: {is_best: {_eq: true}}) {
      edited_positions edited_cds_positions
      post_codon { nt1 nt2 nt3 }
    }
  }
}
"""

INSERT_GUIDES = """
mutation($objects: [guides_insert_input!]!) {
  insert_guides(objects: $objects,
    on_conflict: {constraint: guides_sequence_type_key, update_columns: [hits_85, hits_90, hits_95, hits_100]}
  ) { affected_rows returning { id sequence type } }
}
"""
INSERT_BYSTANDERS = """
mutation($objects: [bystanders_insert_input!]!) {
  insert_bystanders(objects: $objects, on_conflict: {constraint: bystanders_pkey, update_columns: []}) { affected_rows }
}
"""


def _prepare_guide(g):
    return {"sequence": g["sequence"], "type": g["type"],
            "hits_85": g["hits_85"], "hits_90": g["hits_90"], "hits_95": g["hits_95"], "hits_100": g["hits_100"],
            "variants_guides": {"data": [{"variant_id": g["variant_id"], "variant_idx": g["variant_idx"],
                                          "edit_type": g.get("edit_type"), "rescue_sift": None}],
                                "on_conflict": {"constraint": "variants_guides_pkey",
                                                "update_columns": ["variant_idx", "edit_type", "rescue_sift"]}}}


def _fetch_variants(gql, limit):
    rows, offset, page = [], 0, 1000
    while True:
        ps = page if limit is None else min(page, limit - len(rows))
        if ps <= 0:
            break
        batch = gql(GET_IMPROVE_GUIDE_VARIANTS, {"limit": ps, "offset": offset})["variants"]
        rows.extend(batch)
        if len(batch) < ps:
            break
        offset += ps
    return rows


def populate_improve_guides(limit=None, dry_run=False, forks=8):
    """Build + insert Improve guides + bystanders (online CADD/SIFT). Reads paths from env."""
    global genome_fasta, cds_fasta, needs_chr_prefix, BLAT_PATH, BLAT_WORKERS, BLAT_DIR, VEP_SERVER
    gql = hp.gql
    GENOME = _resolve(os.environ["GENOME_REFERENCE_HG19"])
    CDS_REFERENCE = _resolve(os.environ["CDS_REFERENCE"])
    MANE_CDS_REFERENCE = _resolve(os.environ["MANE_CDS_REFERENCE"])
    BLAT_PATH = os.environ["BLAT_PATH"]
    BLAT_WORKERS = int(os.environ["BLAT_WORKERS"])
    VEP_SERVER = os.environ["GRCH37_VEP_SERVER"]
    OUTPUT_DIR = _resolve(os.environ["OUTPUT_DIR"])
    BLAT_DIR = OUTPUT_DIR / "blat_results"
    BLAT_DIR.mkdir(parents=True, exist_ok=True)
    gtf = _resolve(os.environ.get("NEIGHBOR_DRILL_GTF", "Resources/Homo_sapiens.GRCh37.87.gtf.gz"))

    genome_fasta = Fasta(str(GENOME), as_raw=True)
    cds_fasta = Fasta(str(CDS_REFERENCE), as_raw=True)
    needs_chr_prefix = any(k.startswith("chr") for k in genome_fasta.keys())

    codon_lookup = {ov.to_dna(f"{c['nt1']}{c['nt2']}{c['nt3']}"): c["id"]
                    for c in gql("query { codons { id nt1 nt2 nt3 } }", {})["codons"]}
    variants = _fetch_variants(gql, limit)
    log.info(f"[improve_guides] variants: {len(variants)}")
    if not variants:
        return
    vmap = {v["id"]: v for v in variants}
    wanted = {v["variants_features"][0]["feature"]["identifier"]
              for v in variants if v.get("variants_features")}
    tx = ov.load_cds_intervals(gtf, wanted)

    results = [process_variant_cpu(v) for v in variants]

    # ── per-variant codon geometry: intended-edit positions + variant-codon positions ──
    # geom[vid] = {edited_g, edited_c, codon_g(list of 3), codon_cds_start, post_dna, tid_ver, strand}
    geom = {}
    for v in variants:
        vid = v["id"]
        best = (v.get("best") or [None])[0]
        feats = v.get("variants_features") or []
        feat = (feats[0].get("feature") if feats else None) or {}
        tid = feat.get("identifier")
        tid_ver = f"{tid}.{feat['version_number']}" if (tid and feat.get("version_number")) else tid
        t = tx.get(tid)
        pic = v.get("position_in_codon")
        cds_pos = v.get("cds_position")
        if not (best and t and pic and cds_pos):
            geom[vid] = None
            continue
        codon_cds_start = int(cds_pos) - (int(pic) - 1)
        codon_c = [codon_cds_start + i for i in range(3)]
        codon_g = [ov.cds_to_genomic(t, p) for p in codon_c]
        edited_c = set(best.get("edited_cds_positions") or [])
        edited_g = {ov.cds_to_genomic(t, p) for p in edited_c}
        post_dna = ov.to_dna(f"{best['post_codon']['nt1']}{best['post_codon']['nt2']}{best['post_codon']['nt3']}") \
            if best.get("post_codon") else None
        geom[vid] = {"edited_g": edited_g, "edited_c": edited_c, "codon_g": codon_g, "codon_c": codon_c,
                     "codon_cds_start": codon_cds_start, "post_dna": post_dna, "tid_ver": tid_ver,
                     "strand": (v["gene"] or {}).get("strand") or "+"}

    # ── build bystander HGVS (lone A→G for all; codon-delins for same-codon), collect ──
    # each retained bystander gets: kind ('p'/'r'), _lone (hgvs), _same (bool), _delins (hgvs|None), _gpos/_cpos
    all_hgvs = []
    for r in results:
        vid = r["variant_id"]
        g = geom.get(vid)
        if g is None:
            r["pre_bystanders"] = []; r["rna_bystanders"] = []
            continue
        chrom_vep = r["chrom"][3:] if r["chrom"].startswith("chr") else r["chrom"]
        allele = "T>C" if r["strand"] == "-" else "A>G"
        # pre (genomic)
        kept = []
        for b in r["pre_bystanders"]:
            gp = b["position"]
            if gp in g["edited_g"]:
                continue                                   # intended edit → not a bystander
            b["_lone"] = f"{chrom_vep}:g.{gp}{allele}"
            b["_gpos"] = gp
            b["_same"] = gp in g["codon_g"]
            b["_delins"] = _codon_delins(g, codon_idx=(g["codon_g"].index(gp) if b["_same"] else None))
            kept.append(b); all_hgvs.append(b["_lone"])
            if b["_delins"]:
                all_hgvs.append(b["_delins"])
        r["pre_bystanders"] = kept
        # mature (cds)
        kept = []
        for b in r["rna_bystanders"]:
            cp = b["position"]
            if cp in g["edited_c"]:
                continue
            b["_lone"] = f"{g['tid_ver']}:c.{cp}A>G"
            b["_cpos"] = cp
            b["_same"] = cp in g["codon_c"]
            b["_delins"] = _codon_delins(g, codon_idx=(cp - g["codon_cds_start"] if b["_same"] else None))
            kept.append(b); all_hgvs.append(b["_lone"])
            if b["_delins"]:
                all_hgvs.append(b["_delins"])
        r["rna_bystanders"] = kept

    n_pre = sum(len(r["pre_bystanders"]) for r in results)
    n_rna = sum(len(r["rna_bystanders"]) for r in results)
    log.info(f"[improve_guides] bystanders: pre={n_pre} mature={n_rna} | VEP HGVS queries: {len(set(all_hgvs))}")
    vep = _vep_hgvs(all_hgvs) if all_hgvs else {}

    def bystander_obj(r, kind, b):
        g = geom[r["variant_id"]]
        tid = r["transcript_id"]
        lone_tc = _pick_tc(vep.get(b["_lone"]), tid)
        score_tc = _pick_tc(vep.get(b["_delins"]), tid) if b.get("_delins") else lone_tc
        codons = (score_tc.get("codons") or "").upper()
        ref_c = codons.split("/")[0] if "/" in codons else None
        alt_c = codons.split("/")[1] if "/" in codons else None
        cons = score_tc.get("consequence_terms") or []
        return {"in_cds": kind == "r",
                "genomic_position": b["_gpos"] if kind == "p" else None,
                "cds_position": b["position"] if kind == "r" else None,
                "has_u_before": b["has_u_before"], "has_u_after": b["has_u_after"],
                "cadd_raw": lone_tc.get("cadd_raw"), "cadd_phred": lone_tc.get("cadd_phred"),
                "sift_score": score_tc.get("sift_score"), "sift_prediction": score_tc.get("sift_prediction"),
                "ref_codon_id": codon_lookup.get(ref_c) if ref_c else None,
                "alt_codon_id": codon_lookup.get(alt_c) if alt_c else None,
                "consequence_names": hp.consequence_terms_to_regions(cons)}

    all_guides = []
    for r in results:
        if r["dna_guide_obj"]:
            r["dna_guide_obj"]["bystanders"] = [bystander_obj(r, "p", b) for b in r["pre_bystanders"]]
            all_guides.append(r["dna_guide_obj"])
        if r["rna_guide_obj"]:
            r["rna_guide_obj"]["bystanders"] = [bystander_obj(r, "r", b) for b in r["rna_bystanders"]]
            all_guides.append(r["rna_guide_obj"])

    # ── BLAT off-target (guides are around the variant) ──
    pre = [g for g in all_guides if g["type"] == "pre_mRNA" and g["sequence"]]
    mat = [g for g in all_guides if g["type"] == "mature_mRNA" and g["sequence"]]
    pre_hits = blat_all_guides(pre, GENOME, "pre_mRNA")
    mat_hits = blat_all_guides(mat, MANE_CDS_REFERENCE, "mature_mRNA")
    for g in all_guides:
        h = (pre_hits if g["type"] == "pre_mRNA" else mat_hits).get(str(g["variant_id"]))
        for thr in IDENTITY_THRESHOLDS:
            g[f"hits_{thr}"] = (h or {}).get(f"hits_{thr}", 0)

    n_by = sum(len(g["bystanders"]) for g in all_guides)
    n_cadd = sum(1 for g in all_guides for b in g["bystanders"] if b["cadd_phred"] is not None)
    n_same = sum(1 for r in results for b in (r["pre_bystanders"] + r["rna_bystanders"]) if b.get("_delins"))
    log.info(f"[improve_guides] guides: {len(all_guides)} | bystanders: {n_by} (cadd non-null: {n_cadd}) "
             f"| same-codon (codon-aware SIFT): {n_same}")
    if dry_run:
        log.info("[improve_guides] --dry-run: not writing.")
        return

    guide_lookup = {(g["sequence"], g["type"]): g for g in all_guides}
    bystander_objs = []
    for i in tqdm(range(0, len(all_guides), 100), desc="Insert guides"):
        res = gql(INSERT_GUIDES, {"objects": [_prepare_guide(g) for g in all_guides[i:i + 100]]})
        for ret in res["insert_guides"]["returning"]:
            g = guide_lookup.get((ret["sequence"], ret["type"]))
            if not g:
                continue
            for b in g["bystanders"]:
                obj = {"guide_id": ret["id"], "variant_id": g["variant_id"], "in_cds": b["in_cds"],
                       "genomic_position": b["genomic_position"], "cds_position": b["cds_position"],
                       "has_u_before": b["has_u_before"], "has_u_after": b["has_u_after"],
                       "cadd_raw": b["cadd_raw"], "cadd_phred": b["cadd_phred"],
                       "sift_score": b["sift_score"], "sift_prediction": b["sift_prediction"],
                       "ref_codon_id": b["ref_codon_id"], "alt_codon_id": b["alt_codon_id"]}
                if b["consequence_names"]:
                    obj["bystanders_consequences"] = {"data": [
                        {"consequence": {"data": {"name": n}, "on_conflict": {"constraint": "consequences_name_key", "update_columns": ["name"]}}}
                        for n in b["consequence_names"]],
                        "on_conflict": {"constraint": "bystanders_consequences_pkey", "update_columns": []}}
                bystander_objs.append(obj)

    tot_by = 0
    for i in tqdm(range(0, len(bystander_objs), 500), desc="Insert bystanders"):
        tot_by += gql(INSERT_BYSTANDERS, {"objects": bystander_objs[i:i + 500]})["insert_bystanders"]["affected_rows"]
    log.info(f"[improve_guides] ✅ inserted {len(all_guides)} guides and {tot_by} bystander rows (affected).")


def _codon_delins(g, codon_idx):
    """Coding delins for a same-codon bystander: therapy codon (post_codon) + the bystander A→G."""
    if codon_idx is None or g.get("post_dna") is None:
        return None
    chars = list(g["post_dna"])
    if chars[codon_idx] != "A":            # bystander must be an unedited A in the therapy codon
        return None
    chars[codon_idx] = "G"
    full_alt = "".join(chars)
    s = g["codon_cds_start"]
    return f"{g['tid_ver']}:c.{s}_{s + 2}delins{full_alt}"
