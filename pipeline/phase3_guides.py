
# --- imports this stage needs -------------------------------------------------
# Previously inherited from the runner's scope via exec(). A stage that cannot say
# what it imports cannot be read, tested, or run on its own, so each now imports
# for itself; only genuine shared state still arrives from the runner.
import os
import json
import math
import time
import logging
import subprocess
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import Tuple
from pyfaidx import Fasta
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

BLAT_PATH = os.environ["BLAT_PATH"]
CDS_REFERENCE = os.environ["CDS_REFERENCE"]
MANE_CDS_REFERENCE = os.environ["MANE_CDS_REFERENCE"]
GENOME_REFERENCE_HG19 = os.environ["GENOME_REFERENCE_HG19"]
GRCH37_VEP_SERVER = os.environ["GRCH37_VEP_SERVER"]
GRCH38_VEP_SERVER = os.environ["GRCH38_VEP_SERVER"]
VEP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
VEP_REGION_EXT = os.environ["VEP_REGION_EXT"]
VEP_HGVS_EXT = os.environ["VEP_HGVS_EXT"]
BLAT_WORKERS = int(os.environ["BLAT_WORKERS"])

logging.info(f"✅ BLAT configuration loaded:")
logging.info(f"   BLAT: {BLAT_PATH}")
logging.info(f"   MANE SELECT CDS: {MANE_CDS_REFERENCE}")
logging.info(f"   Genome (hg19/GRCh37): {GENOME_REFERENCE_HG19}")
logging.info(f"   Workers: {BLAT_WORKERS}")

# Load FASTA files
cds_fasta = Fasta(CDS_REFERENCE, as_raw=True)
genome_fasta_hg19 = Fasta(GENOME_REFERENCE_HG19, as_raw=True)
genome_fasta_keys = set(genome_fasta_hg19.keys())
needs_chr_prefix = any(k.startswith("chr") for k in genome_fasta_keys)
logging.info(f"✅ FASTA references loaded (needs_chr_prefix={needs_chr_prefix})")


# Step 1: Query Fix variants (G>A)
GET_FIX_VARIANTS = """
query GetFixVariants {
  variants(
    where: {
      _or: [
        { 
          _and: [
            { ref: {_eq: "G"}}, 
            {alt: {_eq: "A"}},
            {_or: [
              {gene: {strand: {_eq: "+"}}},
              {gene: {strand: {_is_null: true}}}
            ]}
          ] 
        }
        { 
          _and: [
            { ref: {_eq: "C"}}, 
            {alt: {_eq: "T"}}, 
            {gene: {strand: {_eq: "-"}}}
          ] 
        }
      ]
      class: {_eq: "SNV"}
    }
  ) {
    id
    ref
    alt
    cds_position
    position_in_codon
    gene_ensg
    gene { strand }
    coordinate { chr start }
    variants_features(where: {feature: {type: {_eq: "Transcript"}}}) {
      feature { identifier version_number }
    }
    ref_codon { id nt1 nt2 nt3 amino_acid { letter } }
    alt_codon { id nt1 nt2 nt3 amino_acid { letter } }
    variants_quantitative_scores(where: {quantitative_score: {name: {_eq: "CADD_phred"}}}) {
      quantitative_score { value }
    }
  }
}
"""

logging.info("Querying Fix (G>A) variants...")
fix_variants = gql(GET_FIX_VARIANTS, {})["variants"]
logging.info(f"✅ Found {len(fix_variants)} Fix variants")



# Step 1: Query Rescue variants (StopGained, excluding G>A)
GET_RESCUE_VARIANTS = """
query GetRescueVariants {
  variants(
    where: {
      variants_consequences: {
        consequence: { name: {_eq: "StopGained"} }
      }
      _not: {
        _or: [
          { 
            _and: [
              { ref: {_eq: "G"}}, 
              {alt: {_eq: "A"}},
              {_or: [
                {gene: {strand: {_eq: "+"}}},
                {gene: {strand: {_is_null: true}}}
              ]}
            ] 
          }
          { 
            _and: [
              { ref: {_eq: "C"}}, 
              {alt: {_eq: "T"}}, 
              {gene: {strand: {_eq: "-"}}}
            ] 
          }
        ]
      }
      class: {_eq: "SNV"}
    }
  ) {
    id
    ref
    alt
    cdna_position
    cds_position
    protein_position
    position_in_codon
    gene_ensg
    gene { strand }
    coordinate { chr start }
    variants_features(where: {feature: {type: {_eq: "Transcript"}}}) {
      feature { 
        identifier 
        version_number
      }
    }
    ref_codon { id nt1 nt2 nt3 amino_acid { letter } }
    alt_codon { id nt1 nt2 nt3 }
    variants_quantitative_scores(where: {quantitative_score: {name: {_eq: "CADD_phred"}}}) {
      quantitative_score { value }
    }
  }
}
"""

logging.info("Querying Rescue (StopGained) variants (excluding G>A)...")
rescue_variants = gql(GET_RESCUE_VARIANTS, {})["variants"]
logging.info(f"✅ Found {len(rescue_variants)} Rescue variants")



def calculate_psl_percent_identity(matches, misMatches, repMatches, qNumInsert, tNumInsert, 
                                   qStart, qEnd, tStart, tEnd, sizeMul=1):
    """
    Calculate percent identity exactly as UCSC BLAT web interface does.
    Based on pslCalcMilliBad from UCSC's `https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/utils/pslScore/pslScore.pl` file.
    """

    qAliSize = sizeMul * (qEnd - qStart)
    tAliSize = tEnd - tStart
    aliSize = min(qAliSize, tAliSize)
    
    if aliSize <= 0:
        return 0.0
    
    sizeDif = abs(qAliSize - tAliSize)
    insertFactor = qNumInsert + tNumInsert
    
    total = sizeMul * (matches + repMatches + misMatches)
    if total == 0:
        return 0.0
    
    # UCSC's log-scaled penalty for size differences
    roundAwayFromZero = 3 * math.log(1 + sizeDif)
    roundAwayFromZero = int(roundAwayFromZero + 0.5 if roundAwayFromZero >= 0 else roundAwayFromZero - 0.5)
    
    milliBad = (1000 * (misMatches * sizeMul + insertFactor + roundAwayFromZero)) / total
    
    return 100.0 - (milliBad * 0.1)



def _extract_guide(fasta_handle: Fasta, key: str, pos: int, strand: str = '+', n: int = 20) -> Tuple[str, int]:
    """
    Extract a guide sequence of length 2n+1 centered around the given position from the FASTA handle.
    Returns:
    - guide_seq: The extracted guide sequence (reverse complemented if strand is '-')
    - variant_idx: The 0-based index of the variant position within the guide sequence
    """
    full_seq = fasta_handle[key]
    seq_len = len(full_seq)
    idx = pos - 1  # Convert to 0-based index

    if idx < 0 or idx >= seq_len:
        raise ValueError(f"Position {pos} out of bounds for sequence {key} of length {seq_len}")
    
    start = max(0, idx - n)
    end = min(seq_len, idx + n + 1)
    guide_seq = str(full_seq[start:end])
    variant_idx = idx - start

    if strand == '-':
        guide_seq = str(Seq(guide_seq).reverse_complement())
        variant_idx = len(guide_seq) - 1 - variant_idx  # Mirror the index after reverse complement

    guide_seq = guide_seq.upper()

    if len(guide_seq) < 2 * n + 1 and pos - n > 0 and pos + n <= seq_len:
        raise ValueError(f"Extracted guide sequence is too short: expected {2*n+1} bases, got {len(guide_seq)} for {key}:{pos} ({strand})")
    
    return guide_seq, variant_idx


GET_CODONS = """
query GetCodons {
  codons {
    id
    nt1
    nt2
    nt3
    amino_acid { letter short_name }
  }
}
"""

codons_data = gql(GET_CODONS, {})["codons"]
codons_lookup = {f"{c['nt1']}{c['nt2']}{c['nt3']}": c for c in codons_data}
logging.info(f"✅ Loaded {len(codons_lookup)} codons")
logging.info(codons_lookup)



def _find_bystanders(guide: str, variant_idx: int, position: int, strand: str = '+', window_size: int = 10):
    """
    Find adenosines within ±window_size of variant position in RNA context.
    
    Args:
        guide: guide sequence (already reverse-complemented for '-' strand)
        variant_idx: Index of the variant in the guide (0-based)
        position: position of the variant (genomic for pre-mRNA, CDS for mature mRNA)
        strand: '+' or '-' — determines how offsets map to positions
        window_size: Window size around variant (default: 10)
    
    Returns:
        List of bystander dicts with all required schema fields
    """
    try:
        position = int(position)
    except (ValueError, TypeError):
        raise ValueError(f"Invalid position value: {position}")
    
    if guide is None:
        raise ValueError("Guide sequence must be provided")
    
    if variant_idx is None:
        raise ValueError("Variant index must be provided")
    
    
    context = guide.upper().replace('T', 'U')
    context_len = len(context)
    
    if variant_idx < 0 or variant_idx >= context_len:
        raise ValueError(f"Variant index {variant_idx} out of bounds for guide of length {context_len}")
    
    bystanders = []
    start_idx = max(0, variant_idx - window_size)
    end_idx = min(context_len, variant_idx + window_size + 1)
    
    for i in range(start_idx, end_idx):
        if i == variant_idx:  # Skip variant position
            continue
        
        if context[i] == 'A':
            # Calculate absolute position accounting for strand
            # '+' strand or CDS: moving right in guide = increasing position
            # '-' strand (genomic): moving right in guide (RC'd) = decreasing genomic position
            offset = i - variant_idx
            if strand == '-':
                offset = -offset
            bystander_pos = position + offset
            
            # Check for adjacent U (context already has U not T)
            has_U_before = (i > 0 and context[i - 1] == 'U')
            has_U_after = (i < context_len - 1 and context[i + 1] == 'U')
            
            bystander_dict = {
                "has_u_before": has_U_before,
                "has_u_after": has_U_after,
                "position": bystander_pos,
            }
            
            bystanders.append(bystander_dict)
    
    return bystanders




# graphql mutation for nested guides inserting 
INSERT_GUIDES = """
mutation InsertGuides($objects: [guide_insert_input!]!) {
  insert_guide(objects: $objects) {
    affected_rows
    returning {
      id
      sequence
      variant_id
      variant { 
        id
        ref
        alt
        coordinate { chr start }
        gene { strand }
      }
      bystanders {
        id
        position
        has_adjacent_u
      }
    }
  }
}
"""




# ── Cell: Helper functions ──

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
TGG = "TGG"


def _build_rescue_sift_query(variant: dict) -> tuple[int, str, str] | None:
    """
    Build a VEP coding HGVS query for ref_codon → TGG using the variant's ENST.
    Uses single-nt substitution when only 1 position differs, delins otherwise.

    Returns:
      (variant_id, hgvs_coding, transcript_id) — needs VEP lookup
      None                                     — ref is already TGG or missing data
    """
    ref_codon_obj = variant.get("ref_codon")
    if not ref_codon_obj:
        return None
    ref_codon = f"{ref_codon_obj['nt1']}{ref_codon_obj['nt2']}{ref_codon_obj['nt3']}"
    if ref_codon == TGG:
        return None

    cds_position = variant.get("cds_position")
    position_in_codon = variant.get("position_in_codon")
    if cds_position is None or position_in_codon is None:
        return None

    variants_features = variant.get("variants_features", [])
    if not variants_features:
        return None
    feature = variants_features[0].get("feature")
    if not feature:
        return None

    transcript_id = feature.get("identifier")
    version = feature.get("version_number")
    if not transcript_id:
        return None

    enst_ver = f"{transcript_id}.{version}" if version else transcript_id

    cds_pos = int(cds_position)
    pic = int(position_in_codon)  # 1-based
    codon_cds_start = cds_pos - (pic - 1)

    # Find which codon positions differ from TGG
    target = "TGG"
    diffs = [(i, ref_codon[i], target[i]) for i in range(3) if ref_codon[i] != target[i]]

    if len(diffs) == 1:
        # Single nucleotide substitution → SIFT will always work
        i, ref_nt, alt_nt = diffs[0]
        pos = codon_cds_start + i
        hgvs = f"{enst_ver}:c.{pos}{ref_nt}>{alt_nt}"
    else:
        # Multi-nucleotide: use delins
        codon_cds_end = codon_cds_start + 2
        hgvs = f"{enst_ver}:c.{codon_cds_start}_{codon_cds_end}delinsTGG"

    return (variant["id"], hgvs, transcript_id)

def _extract_sift_from_vep(entry: dict | None) -> float | None:
    """Extract sift_score from a VEP region result entry."""
    if entry is None:
        return None
    tcs = entry.get("transcript_consequences", [])
    if tcs:
        return tcs[0].get("sift_score")
    return None


def _vep_region_batch(region_list: list[str], batch_size: int = 200, verbose: bool = False, max_workers: int = 8) -> list[dict | None]:
    """POST batches of VEP region strings, return results in same order."""
    num_batches = math.ceil(len(region_list) / batch_size)
    all_results: list[dict | None] = [None] * len(region_list)

    def _send_batch(start: int) -> tuple[int, list[dict | None]]:
        chunk = region_list[start : start + batch_size]
        payload = json.dumps({"variants": chunk})
        for attempt in range(5):
            r = requests.post(
                GRCH37_VEP_SERVER + "/vep/human/region?CADD=1&pick=1",
                headers=VEP_HEADERS, data=payload, timeout=120,
            )
            if r.status_code == 200:
                break
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 1)))
                continue
            if r.status_code in (502, 503, 504):
                time.sleep(0.5 * (2 ** attempt))
                continue
            raise RuntimeError(f"VEP region {r.status_code}: {r.text[:800]}")
        else:
            raise RuntimeError(f"VEP region failed after retries (batch at {start})")

        vep_responses = r.json()
        batch_results = []
        for entry in vep_responses:
            batch_results.append(entry)
        while len(batch_results) < len(chunk):
            batch_results.append(None)
        return start, batch_results

    batch_starts = list(range(0, len(region_list), batch_size))

    with ThreadPoolExecutor(max_workers=min(max_workers, num_batches)) as pool:
        futures = {pool.submit(_send_batch, s): s for s in batch_starts}
        completed = as_completed(futures)
        if verbose:
            completed = tqdm(completed, total=num_batches, desc="VEP region batches")
        for future in completed:
            start, batch_results = future.result()
            for i, res in enumerate(batch_results):
                all_results[start + i] = res

    return all_results


def _make_guide_obj(guide_type: str, variant_id: int, sequence: str, variant_idx: int, edit_type: str = None, rescue_sift: float | None = None) -> dict:
    return {
        "hits_85": None, "hits_90": None, "hits_95": None, "hits_100": None,
        "sequence": sequence,
        "variant_id": variant_id,
        "variant_idx": variant_idx,
        "bystanders": [],
        "type": guide_type,
        "edit_type": edit_type,
        "rescue_sift": rescue_sift,
    }


def _vep_batch(hgvs_list: list[str], batch_size: int = 200, verbose: bool = False, max_workers: int = 8) -> list[dict | None]:
    """POST up to `batch_size` HGVS notations per request, with parallel batch submission."""
    num_batches = math.ceil(len(hgvs_list) / batch_size)
    all_results: list[dict | None] = [None] * len(hgvs_list)

    def _send_batch(start: int) -> tuple[int, list[dict | None]]:
        chunk = hgvs_list[start : start + batch_size]
        payload = json.dumps({"hgvs_notations": chunk})
        for attempt in range(5):
            r = requests.post(
                GRCH37_VEP_SERVER + VEP_HGVS_EXT,
                headers=VEP_HEADERS, data=payload, timeout=120,
            )
            if r.status_code == 200:
                break
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 1)))
                continue
            if r.status_code in (502, 503, 504):
                time.sleep(0.5 * (2 ** attempt))
                continue
            raise RuntimeError(f"VEP {r.status_code}: {r.text[:800]}")
        else:
            raise RuntimeError(f"VEP failed after retries (batch at {start})")

        results_by_input = {e["input"]: e for e in r.json()}
        return start, [results_by_input.get(hgvs) for hgvs in chunk]

    batch_starts = list(range(0, len(hgvs_list), batch_size))

    with ThreadPoolExecutor(max_workers=min(max_workers, num_batches)) as pool:
        futures = {pool.submit(_send_batch, s): s for s in batch_starts}
        completed = as_completed(futures)
        if verbose:
            completed = tqdm(completed, total=num_batches, desc="VEP batches")
        for future in completed:
            start, batch_results = future.result()
            for i, res in enumerate(batch_results):
                all_results[start + i] = res

    return all_results


def _pick_tc(entry: dict | None, transcript_id: str | None = None) -> dict:
    if entry is None:
        return {}
    tcs = entry.get("transcript_consequences", [])
    if transcript_id:
        return next((t for t in tcs if t.get("transcript_id") == transcript_id), {})
    return tcs[0] if tcs else {}


def _build_bystander_obj(bystander: dict, tc: dict, in_cds: bool) -> dict:
    am = tc.get("alphamissense") or {}
    codons = tc.get("codons")
    ref_codon = codons.split("/")[0].upper() if codons else None
    alt_codon = codons.split("/")[1].upper() if codons else None
    return {
        "in_cds": in_cds,
        "genomic_position": bystander["position"] if not in_cds else None,
        "cds_position": bystander["position"] if in_cds else None,
        "has_u_before": bystander["has_u_before"],
        "has_u_after": bystander["has_u_after"],
        "am_pathogenicity": am.get("am_pathogenicity"),
        "am_class": am.get("am_class"),
        "sift_score": tc.get("sift_score"),
        "sift_prediction": tc.get("sift_prediction"),
        "cadd_raw": tc.get("cadd_raw"),
        "cadd_phred": tc.get("cadd_phred"),
        "ref_codon_id": codons_lookup.get(ref_codon, {}).get("id") if ref_codon else None,
        "alt_codon_id": codons_lookup.get(alt_codon, {}).get("id") if alt_codon else None,
        # Region names for bystanders_consequences (mirrors variants_consequences).
        "consequence_names": hp.consequence_terms_to_regions(tc.get("consequence_terms")),
    }

def _process_variant_cpu(variant: dict, edit_type: str = None) -> dict | None:
    """
    CPU-only work per variant: extract guides + find bystanders + build HGVS strings.
    Returns a dict with everything needed for the VEP stage, or None on failure.
    No network calls here — safe for threads/processes.
    """
    vid = variant["id"]
    strand = (variant["gene"] or {}).get("strand") or "+"
    chrom = str(variant["coordinate"]["chr"])
    if needs_chr_prefix and not chrom.startswith("chr"):
        chrom = "chr" + chrom
    gpos = variant["coordinate"]["start"]

    result = {"variant_id": vid, "errors": [], "edit_type": edit_type}

    # ── pre-mRNA guide ──
    try:
        dna_guide, dna_vi = _extract_guide(genome_fasta_hg19, chrom, gpos, strand, 20)
        result["dna_guide_obj"] = _make_guide_obj("pre_mRNA", vid, dna_guide, dna_vi + 1, edit_type=edit_type)
    except Exception as e:
        result["errors"].append((vid, "pre_mRNA", str(e)))
        return result  # can't continue without pre-mRNA guide

    # ── mature mRNA guide ──
    # Initialize as None - will be populated if possible
    result["rna_guide_obj"] = None
    result["transcript_id"] = None
    result["tid_ver"] = None
    
    # check if the variant has transcript info and a valid CDS position
    if variant.get("variants_features") and variant["variants_features"][0].get("feature") and variant.get("cds_position"):
        try:
            transcript_id = variant["variants_features"][0]["feature"]["identifier"]
            tid_ver = f'{transcript_id}.{variant["variants_features"][0]["feature"]["version_number"]}'
            cds_pos = int(variant["cds_position"])
            rna_guide, rna_vi = _extract_guide(cds_fasta, tid_ver, cds_pos)
            result["rna_guide_obj"] = _make_guide_obj("mature_mRNA", vid, rna_guide, rna_vi + 1, edit_type=edit_type)
            result["transcript_id"] = transcript_id
            result["tid_ver"] = tid_ver
        except Exception as e:
            # Log the error but continue - pre-mRNA guide is still valid
            result["errors"].append((vid, "mature_mRNA", str(e)))
            result["rna_guide_obj"] = None

    # ── bystanders (CPU) ──
    allele_string = "T>C" if strand == "-" else "A>G"

    # Process pre-mRNA bystanders
    try:
        pre_bystanders = _find_bystanders(dna_guide, dna_vi, gpos, strand)
        pre_hgvs = [f"{chrom}:g.{b['position']}{allele_string}" for b in pre_bystanders]
        result["pre_bystanders"] = pre_bystanders
        result["pre_hgvs"] = pre_hgvs
    except Exception as e:
        result["errors"].append((vid, "pre_mRNA_bystanders", str(e)))
        result["pre_bystanders"] = []
        result["pre_hgvs"] = []

    # Only process RNA bystanders if we have an RNA guide
    if result["rna_guide_obj"] is not None:
        try:
            rna_guide = result["rna_guide_obj"]["sequence"]
            cds_pos = int(variant["cds_position"])
            
            rna_bystanders = _find_bystanders(rna_guide, rna_vi, cds_pos)
            rna_hgvs = [f"{tid_ver}:c.{b['position']}A>G" for b in rna_bystanders]
            
            result["rna_bystanders"] = rna_bystanders
            result["rna_hgvs"] = rna_hgvs
        except Exception as e:
            result["errors"].append((vid, "mature_mRNA_bystanders", str(e)))
            result["rna_bystanders"] = []
            result["rna_hgvs"] = []
    else:
        result["rna_bystanders"] = []
        result["rna_hgvs"] = []

    return result



# ── Cell: Main pipeline (parallel CPU + batched VEP) ──

VERBOSE = True  # Set to False to suppress all progress bars and prints

all_variants = fix_variants + rescue_variants
fix_variant_ids = {v["id"] for v in fix_variants}
rescue_variant_ids = {v["id"] for v in rescue_variants}
rescue_variant_lookup = {v["id"]: v for v in rescue_variants}

failed_variants = []
all_guide_objects = []

# ═══════════════════════════════════════════════════════
# Phase 1: Parallel CPU work (guide extraction + bystanders)
# ═══════════════════════════════════════════════════════
if VERBOSE:
    logging.info(f"Phase 1: Extracting guides & bystanders for {len(all_variants)} variants...")

with ThreadPoolExecutor(max_workers=min(32, os.cpu_count() or 4)) as pool:
    futures = {}
    for i, v in enumerate(all_variants):
        et = "Fix" if v["id"] in fix_variant_ids else "Rescue"
        futures[pool.submit(_process_variant_cpu, v, et)] = i
    cpu_results = [None] * len(all_variants)
    completed = as_completed(futures)
    if VERBOSE:
        completed = tqdm(completed, total=len(futures), desc="CPU extraction")
    for future in completed:
        idx = futures[future]
        try:
            cpu_results[idx] = future.result()
        except Exception as e:
            cpu_results[idx] = None
            failed_variants.append((all_variants[idx]["id"], "cpu_exception", str(e)))

# Collect errors from CPU phase & filter valid results
valid_results = []
for res in cpu_results:
    if res is None:
        continue
    failed_variants.extend(res.get("errors", []))
    # Must have dna_guide_obj (always required) and rna_guide_obj key (even if None)
    if "dna_guide_obj" in res:
        valid_results.append(res)

if VERBOSE:
    logging.info(f"  ✅ {len(valid_results)} variants ready for VEP, {len(failed_variants)} failures so far")

# ═══════════════════════════════════════════════════════
# Phase 2: Collect ALL HGVS strings across variants → one mega-batch VEP call
# ═══════════════════════════════════════════════════════
if VERBOSE:
    logging.info("Phase 2: Batched VEP annotation...")

global_hgvs = []
bookkeeping = []

for res in valid_results:
    start = len(global_hgvs)
    n_pre = len(res["pre_hgvs"])
    n_rna = len(res["rna_hgvs"])
    global_hgvs.extend(res["pre_hgvs"])
    global_hgvs.extend(res["rna_hgvs"])
    bookkeeping.append((start, n_pre, n_rna))

if VERBOSE:
    logging.info(f"  Total HGVS queries: {len(global_hgvs)}")

if global_hgvs:
    try:
        global_vep_results = _vep_batch(global_hgvs, batch_size=200, verbose=VERBOSE)
    except Exception as e:
        if VERBOSE:
            logging.info(f"  ❌ VEP batch failed: {e}")
        global_vep_results = [None] * len(global_hgvs)
else:
    global_vep_results = []

# ═══════════════════════════════════════════════════════
# Phase 2.5: Rescue SIFT — VEP protein HGVS queries for ref AA→Trp
# ═══════════════════════════════════════════════════════
if VERBOSE:
    logging.info("Phase 2.5: Rescue SIFT (protein-level ref→Trp VEP queries)...")

rescue_sift_queries = []  # list of (variant_id, hgvs_protein, transcript_id)
rescue_sift_lookup = {}   # {variant_id: sift_score}

for vid, v in rescue_variant_lookup.items():
    result = _build_rescue_sift_query(v)
    if isinstance(result, tuple):
        rescue_sift_queries.append(result)
    # else None → missing data, skip

if rescue_sift_queries:
    rescue_vids = [q[0] for q in rescue_sift_queries]
    rescue_hgvs = [q[1] for q in rescue_sift_queries]
    rescue_transcripts = [q[2] for q in rescue_sift_queries]

    if VERBOSE:
        logging.info(f"  Querying VEP for {len(rescue_hgvs)} rescue protein-level changes...")

    try:
        rescue_vep_results = _vep_batch(rescue_hgvs, batch_size=200, verbose=VERBOSE)
        for vid, entry, transcript_id in zip(rescue_vids, rescue_vep_results, rescue_transcripts):
            # Extract SIFT using the specific transcript
            tc = _pick_tc(entry, transcript_id)
            rescue_sift_lookup[vid] = tc.get("sift_score")
    except Exception as e:
        if VERBOSE:
            logging.info(f"  ❌ Rescue VEP batch failed: {e}")

if VERBOSE:
    n_with_sift = sum(1 for v in rescue_sift_lookup.values() if v is not None)
    logging.info(f"  ✅ {len(rescue_sift_queries)} rescue variants queried, {n_with_sift} got SIFT scores")

# ═══════════════════════════════════════════════════════
# Phase 3: Assemble bystander objects from VEP results + attach rescue SIFT
# ═══════════════════════════════════════════════════════
if VERBOSE:
    logging.info("Phase 3: Building guide objects...")

assembly_iter = zip(valid_results, bookkeeping)
if VERBOSE:
    assembly_iter = tqdm(assembly_iter, total=len(valid_results), desc="Assembling guides")

for res, (start, n_pre, n_rna) in assembly_iter:
    vid = res["variant_id"]
    pre_vep = global_vep_results[start : start + n_pre]
    rna_vep = global_vep_results[start + n_pre : start + n_pre + n_rna]

    try:
        res["dna_guide_obj"]["bystanders"] = [
            _build_bystander_obj(b, _pick_tc(entry), in_cds=False)
            for b, entry in zip(res["pre_bystanders"], pre_vep)
        ]
    except Exception as e:
        failed_variants.append((vid, "pre_mRNA_bystanders", str(e)))

    # Only process RNA bystanders if RNA guide exists
    if res["rna_guide_obj"] is not None:
        try:
            res["rna_guide_obj"]["bystanders"] = [
                _build_bystander_obj(b, _pick_tc(entry, res["transcript_id"]), in_cds=True)
                for b, entry in zip(res["rna_bystanders"], rna_vep)
            ]
        except Exception as e:
            failed_variants.append((vid, "mature_mRNA_bystanders", str(e)))

    # ── Attach rescue SIFT (None for Fix, SIFT score for Rescue) ──
    edit_type = res["edit_type"]
    rescue_sift = None

    if edit_type == "Rescue":
        rescue_sift = rescue_sift_lookup.get(vid)

    res["dna_guide_obj"]["rescue_sift"] = rescue_sift
    if res["rna_guide_obj"] is not None:
        res["rna_guide_obj"]["rescue_sift"] = rescue_sift

    # Add guides to insertion list
    all_guide_objects.append(res["dna_guide_obj"])
    if res["rna_guide_obj"] is not None:
        all_guide_objects.append(res["rna_guide_obj"])

if VERBOSE:
    logging.info(f"✅ Done! {len([g for g in all_guide_objects if g['sequence'] is not None])} guide objects, {len(failed_variants)} failures")
    rescue_sifts = [g["rescue_sift"] for g in all_guide_objects if g["edit_type"] == "Rescue" and g["rescue_sift"] is not None]
    logging.info(f"  Rescue SIFT: {len(rescue_sifts)} Rescue guides with non-null SIFT score")
    if failed_variants:
        for vid, step, *msg in failed_variants:
            logging.info(f"  Variant {vid} failed at {step}: {msg[0] if msg else ''}")



# Where are the nulls coming from?
no_ref_codon = sum(1 for v in rescue_variants if not v.get("ref_codon"))
no_pos = sum(1 for v in rescue_variants if v.get("position_in_codon") is None)

built = [_build_rescue_sift_query(v) for v in rescue_variants]
queries_built = sum(1 for b in built if b is not None)
queries_none = sum(1 for b in built if b is None)

logging.info(f"no ref_codon: {no_ref_codon}")
logging.info(f"no position_in_codon: {no_pos}")
logging.info(f"_build_rescue_sift_query returned query: {queries_built}")
logging.info(f"_build_rescue_sift_query returned None: {queries_none}")
logging.info(f"rescue_sift_lookup entries: {len(rescue_sift_lookup)}")
logging.info(f"rescue_sift_lookup non-null: {sum(1 for v in rescue_sift_lookup.values() if v is not None)}")




# Check what VEP actually returned for the failed cases
none_entries = sum(1 for v in rescue_sift_lookup.values() if v is None)
logging.info(f"rescue_sift_lookup None values: {none_entries}")

# Re-check a sample of failed ones
failed_vids = [vid for vid, sift in rescue_sift_lookup.items() if sift is None][:5]
sample_queries = []
for vid in failed_vids:
    v = rescue_variant_lookup[vid]
    result = _build_rescue_sift_query(v)
    if result:
        sample_queries.append(result)
        logging.info(f"  vid={vid}, hgvs={result[1]}, transcript={result[2]}")

# Re-run just these samples through VEP to see what comes back
if sample_queries:
    hgvs_list = [q[1] for q in sample_queries]
    transcripts = [q[2] for q in sample_queries]
    results = _vep_batch(hgvs_list, batch_size=5, verbose=False)
    for q, res, transcript_id in zip(sample_queries, results, transcripts):
        logging.info(f"\n  vid={q[0]}: {q[1]} (transcript={transcript_id})")
        if res is None:
            logging.info(f"    VEP returned: None")
        else:
            tc = _pick_tc(res, transcript_id)
            if tc:
                logging.info(f"    transcript_consequences sift_score={tc.get('sift_score')}, sift_prediction={tc.get('sift_prediction')}")
                logging.info(f"    polyphen_score={tc.get('polyphen_score')}, consequence_terms={tc.get('consequence_terms')}")
            else:
                logging.info(f"    No matching transcript_consequence found")



# ── Cell: Phase 4 — BLAT off-target analysis ──

BLAT_CHUNK_SIZE = 500  # guides per BLAT query file

col_names = [
    'matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert',
    'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd',
    'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts'
]

IDENTITY_THRESHOLDS = [100, 95, 90, 85]


def _run_blat_chunk(chunk_guides: list[dict], reference: str, chunk_id: int, guide_type: str) -> pd.DataFrame:
    """Write a chunk of guides to a temp FASTA, run BLAT, return parsed PSL DataFrame."""
    query_path = BLAT_RESULTS_DIR / f"blat_{guide_type}_chunk{chunk_id}.fa"
    results_path = BLAT_RESULTS_DIR / f"blat_{guide_type}_chunk{chunk_id}.psl"

    # Write query FASTA — use variant_id as the sequence name for later joining
    with open(query_path, 'w') as f:
        for g in chunk_guides:
            f.write(f">{g['variant_id']}\n{g['sequence']}\n")

    blat_cmd = [
        BLAT_PATH, reference, query_path, results_path,
        '-out=psl', '-stepSize=5', '-repMatch=2253',
        '-minScore=20', '-minIdentity=0',
    ]
    subprocess.run(blat_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Parse results
    try:
        df = pd.read_csv(results_path, sep='\t', skiprows=5, names=col_names)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=col_names)

    # Clean up temp files
    os.remove(query_path)
    os.remove(results_path)

    return df


def _blat_all_guides(guides: list[dict], reference: str, guide_type: str, verbose: bool = False) -> dict[int, dict]:
    """
    Run BLAT for a list of guide objects against `reference`.
    Returns {variant_id: {"hits_85": N, "hits_90": N, "hits_95": N, "hits_100": N}}.
    """
    if not guides:
        return {}

    # Split into chunks
    chunks = [guides[i:i + BLAT_CHUNK_SIZE] for i in range(0, len(guides), BLAT_CHUNK_SIZE)]

    # Run BLAT chunks in parallel
    all_dfs = []
    with ProcessPoolExecutor(max_workers=min(BLAT_WORKERS, len(chunks))) as pool:
        futures = {
            pool.submit(_run_blat_chunk, chunk, reference, idx, guide_type): idx
            for idx, chunk in enumerate(chunks)
        }
        completed = as_completed(futures)
        if verbose:
            completed = tqdm(completed, total=len(futures), desc=f"BLAT {guide_type}")
        for future in completed:
            all_dfs.append(future.result())

    if not all_dfs:
        return {}

    blat_df = pd.concat(all_dfs, ignore_index=True)

    # Calculate percent identity
    blat_df['identity'] = blat_df.apply(
        lambda row: calculate_psl_percent_identity(
            row['matches'], row['misMatches'], row['repMatches'],
            row['qNumInsert'], row['tNumInsert'],
            row['qStart'], row['qEnd'], row['tStart'], row['tEnd']
        ), axis=1
    )

    # Convert qName to int (variant_id)
    blat_df['qName'] = blat_df['qName'].astype(str)

    # Count hits per variant_id per threshold
    hit_counts = {}
    for threshold in IDENTITY_THRESHOLDS:
        threshold_df = blat_df[blat_df['identity'] >= threshold]
        counts = threshold_df['qName'].value_counts().to_dict()
        for vid_str, count in counts.items():
            hit_counts.setdefault(vid_str, {})[f"hits_{threshold}"] = max(0, count - 1)  # Subtract 1 to exclude the perfect match to the intended target

    # Fill missing thresholds with 0
    for vid_str in hit_counts:
        for threshold in IDENTITY_THRESHOLDS:
            hit_counts[vid_str].setdefault(f"hits_{threshold}", 0)

    return hit_counts


# ── Split guides by type ──
pre_mrna_guides = [g for g in all_guide_objects if g["type"] == "pre_mRNA" and g["sequence"]]
mature_mrna_guides = [g for g in all_guide_objects if g["type"] == "mature_mRNA" and g["sequence"]]

if VERBOSE:
    logging.info(f"Phase 4: BLAT off-target analysis")
    logging.info(f"  pre_mRNA guides:    {len(pre_mrna_guides)} → {GENOME_REFERENCE_HG19}")
    logging.info(f"  mature_mRNA guides: {len(mature_mrna_guides)} → {MANE_CDS_REFERENCE}")

# ── Run BLAT ──
pre_mrna_hits = _blat_all_guides(pre_mrna_guides, GENOME_REFERENCE_HG19, "pre_mRNA", verbose=VERBOSE)
mature_mrna_hits = _blat_all_guides(mature_mrna_guides, MANE_CDS_REFERENCE, "mature_mRNA", verbose=VERBOSE)


# ── Attach hit counts back to guide objects ──
for g in all_guide_objects:
    vid_str = str(g["variant_id"])
    if g["type"] == "pre_mRNA" and vid_str in pre_mrna_hits:
        g.update(pre_mrna_hits[vid_str])
    elif g["type"] == "mature_mRNA" and vid_str in mature_mrna_hits:
        g.update(mature_mrna_hits[vid_str])
    else:
        # No BLAT hits at all → 0 for all thresholds
        for threshold in IDENTITY_THRESHOLDS:
            g[f"hits_{threshold}"] = 0


if VERBOSE:
    # Summary
    for guide_type, hits_dict in [("pre_mRNA", pre_mrna_hits), ("mature_mRNA", mature_mrna_hits)]:
        n_with_hits = sum(1 for h in hits_dict.values() if h.get("hits_85", 0) > 0)
        logging.info(f"  {guide_type}: {n_with_hits}/{len(hits_dict)} guides have ≥1 hit at 85% identity")
    logging.info(f"✅ BLAT complete. All guide objects updated with hit counts.")




# ── Cell: Phase 5 — Insert guides + bystanders into DB ──

INSERT_GUIDES_MUTATION = """
mutation InsertGuides($objects: [guides_insert_input!]!) {
  insert_guides(
    objects: $objects,
    on_conflict: {
      constraint: guides_sequence_type_key,
      update_columns: [hits_85, hits_90, hits_95, hits_100]
    }
  ) {
    affected_rows
    returning {
      id
      sequence
      type
    }
  }
}
"""

INSERT_BYSTANDERS_MUTATION = """
mutation InsertBystanders($objects: [bystanders_insert_input!]!) {
  insert_bystanders(
    objects: $objects,
    on_conflict: {
      constraint: bystanders_pkey,
      update_columns: []
    }
  ) {
    affected_rows
  }
}
"""

GUIDE_BATCH_SIZE = 100


def _prepare_guide_for_insert(g: dict) -> dict:
    return {
        "sequence": g["sequence"],
        "type": g["type"],
        "hits_85": g["hits_85"],
        "hits_90": g["hits_90"],
        "hits_95": g["hits_95"],
        "hits_100": g["hits_100"],
        "variants_guides": {
            "data": [{"variant_id": g["variant_id"], "variant_idx": g["variant_idx"], "edit_type": g.get("edit_type"), "rescue_sift": g.get("rescue_sift")}],
            "on_conflict": {
                "constraint": "variants_guides_pkey",
                "update_columns": ["variant_idx", "edit_type", "rescue_sift"]
            }
        },
    }


total_guides_inserted = 0
all_bystanders_to_insert = []

batches = [all_guide_objects[i:i + GUIDE_BATCH_SIZE] for i in range(0, len(all_guide_objects), GUIDE_BATCH_SIZE)]

logging.info(f"Phase 5: Inserting {len(all_guide_objects)} guides in {len(batches)} batches...")

# Build a mapping from (sequence, type) → guide object for later bystander linking
guide_lookup = {(g["sequence"], g["type"]): g for g in all_guide_objects}

for batch_num, batch in enumerate(tqdm(batches, desc="Inserting guides")):
    objects = [_prepare_guide_for_insert(g) for g in batch]
    result = gql(INSERT_GUIDES_MUTATION, {"objects": objects})
    total_guides_inserted += result["insert_guides"]["affected_rows"]
    
    # Map returned guide IDs to bystanders
    for returned_guide in result["insert_guides"]["returning"]:
        guide_id = returned_guide["id"]
        sequence = returned_guide["sequence"]
        gtype = returned_guide["type"]
        
        # Find matching guide object
        guide_obj = guide_lookup.get((sequence, gtype))
        if not guide_obj:
            continue
            
        variant_id = guide_obj["variant_id"]
        
        # Prepare bystanders for this guide
        for b in guide_obj.get("bystanders", []):
            bystander_obj = {
                "guide_id": guide_id,
                "variant_id": variant_id,
                "in_cds": b["in_cds"],
                "genomic_position": b.get("genomic_position"),
                "cds_position": b.get("cds_position"),
                "has_u_before": b["has_u_before"],
                "has_u_after": b["has_u_after"],
                "am_pathogenicity": b.get("am_pathogenicity"),
                "am_class": b.get("am_class"),
                "sift_score": b.get("sift_score"),
                "sift_prediction": b.get("sift_prediction"),
                "cadd_raw": b.get("cadd_raw"),
                "cadd_phred": b.get("cadd_phred"),
                "ref_codon_id": b.get("ref_codon_id"),
                "alt_codon_id": b.get("alt_codon_id"),
            }

            # Nested insert for the bystanders_consequences join table — same
            # pattern as variants_consequences in phase2: upsert each consequence
            # by name, then link it to this bystander.
            consequence_names = b.get("consequence_names") or []
            if consequence_names:
                bystander_obj["bystanders_consequences"] = {
                    "data": [
                        {
                            "consequence": {
                                "data": {"name": name},
                                "on_conflict": {
                                    "constraint": "consequences_name_key",
                                    "update_columns": ["name"],
                                },
                            }
                        }
                        for name in consequence_names
                    ],
                    "on_conflict": {
                        "constraint": "bystanders_consequences_pkey",
                        "update_columns": [],
                    },
                }

            all_bystanders_to_insert.append(bystander_obj)

logging.info(f"\n✅ Total guide rows inserted: {total_guides_inserted}")

# Insert bystanders in batches
if all_bystanders_to_insert:
    logging.info(f"\nPhase 6: Inserting {len(all_bystanders_to_insert)} bystanders...")
    bystander_batches = [all_bystanders_to_insert[i:i + 500] for i in range(0, len(all_bystanders_to_insert), 500)]
    total_bystanders = 0
    for batch in tqdm(bystander_batches, desc="Inserting bystanders"):
        result = gql(INSERT_BYSTANDERS_MUTATION, {"objects": batch})
        total_bystanders += result["insert_bystanders"]["affected_rows"]
    logging.info(f"✅ Total bystander rows inserted: {total_bystanders}")



logging.info(f"Finished processing variants{' successfully' if not failed_variants else ''}!")
if failed_variants:
    logging.info(f"However, some variants failed at certain steps:")
    for variant_id, step, *msg in failed_variants:
        error_msg = msg[0] if msg else "Unknown error"
        logging.info(f"  Variant {variant_id} failed at {step}: {error_msg}")



CDS_GROUP = ["StopGained", "Missense", "Splice", "Frameshift", "Synonymous", "StartLost", "StopRetained"]

colors = {
    'cornflowerblue': '#6495ED',
    'orange': '#FFA500',
    'red2': '#EE0000',
    'forestgreen': '#228B22',
    'blue': '#00008B',
    'tan2': '#EE9A49',
    'yellowgreen': '#9ACD32',
    'mediumorchid1': '#E066FF',
    'cadetblue1': '#98F5FF',
    'tan4': '#8B5A2B',
    'black': '#000000',
    'darkseagreen4': '#698B69',
}

cols_consequence = {
  "Missense": colors['orange'],
  "StopGained": colors['red2'],
  "Splice": colors['cornflowerblue'],
  "StartLost": colors['blue'],
  "Synonymous": colors['forestgreen'],
  "Frameshift": colors['tan2'],
  "StopRetained": colors['yellowgreen'],
  "Other": colors['mediumorchid1'],
}


# ── Fetch variants with guide hits and consequences ──

GET_VARIANTS_WITH_HITS = """
query GetVariantsWithHits ($cds_group: [String!]) {
  variants(
    where: {
      class: {_eq: "SNV"}
      _or: [
        {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
        {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
      ]
      variants_features: {
        feature: {biotype: {name: {_eq: "protein_coding"}}}
      }
      gene: {
        genes_quantitative_scores: {
          quantitative_score: {name: {_eq: "SFARI Gene Score"}}
        }
      }
      variants_consequences: {
        consequence: {name: {_in: $cds_group}}
      }
    }
  ) {
    id
    variants_guides {
      guide {
        hits_85
        hits_90
        hits_95
        hits_100
        type
      }
    }
    variants_consequences {
      consequence { name }
    }
  }
}
"""



def plot_hits_consequences(variants_with_hits: list[dict], guide_type: str, hits_threshold: int = 85):


    # ── Parse data ──

    rows = []
    for v in variants_with_hits:
        # Get hits_85 from guide (use first pre_mRNA guide)
        variants_guides = v.get("variants_guides", [])  
        if not variants_guides:
            raise ValueError(f"Variant {v['id']} has no associated guides, but was expected to have at least one.")
        # get the right guide_type
        variant_guide = None
        for vg in variants_guides:
            if vg["guide"]["type"] == guide_type:
                variant_guide = vg["guide"]
                break

        if not variant_guide:
            raise ValueError(f"Variant {v['id']} has no associated guide of type {guide_type}, but was expected to have one.")


        if variant_guide.get(f"hits_{hits_threshold}") is None:  # ← Access through "guide" object
            hits = 0
        else:
            hits = variant_guide[f"hits_{hits_threshold}"]  # ← Access through "guide" object
        
        # Bin: 0, 1, 2, 3, 4, 5+
        hits_bin = str(hits) if hits < 5 else "5+"
        
        # Get consequence (pick the one matching our filter)
        consequence = None
        for vc in v.get("variants_consequences", []):
            consequence = vc["consequence"]["name"]
        
        if consequence:
            rows.append({"hits_bin": hits_bin, "consequence": consequence})

    df_plot = pd.DataFrame(rows)
    logging.info(f"Parsed {len(df_plot)} rows")

    # ── Build stacked bar data ──
    bin_order = ["0", "1", "2", "3", "4", "5+"]
    consequence_order = ["Synonymous", "Missense", "StopGained", "Splice"]

    pivot = df_plot.groupby(["hits_bin", "consequence"]).size().unstack(fill_value=0)
    pivot = pivot.reindex(index=bin_order, columns=consequence_order, fill_value=0)

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(8, 5))

    bottoms = np.zeros(len(bin_order))
    for cons in consequence_order:
        values = pivot[cons].values
        ax.bar(bin_order, values, bottom=bottoms, label=cons, color=cols_consequence[cons], edgecolor='black', linewidth=0.5)
        bottoms += values

    ax.set_xlabel(f"Off-target hits (≥{hits_threshold}% identity)", fontsize=12)
    ax.set_ylabel("Number of variants", fontsize=12)
    ax.set_title("CDS+Splice G>A Variants by Off-Target Hits and Consequence", fontsize=13)
    ax.legend(title="Consequence", frameon=True)

    # Add count labels on top of each bar
    for i, bin_label in enumerate(bin_order):
        total = int(bottoms[i])
        if total > 0:
            ax.text(i, total + max(bottoms) * 0.01, str(total), ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_ylim(0, max(bottoms) * 1.08)
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/plots/off_target_hits_{guide_type}_{hits_threshold}.svg", format="svg", bbox_inches="tight")
    plt.show()

    logging.info(f"\nTotal variants: {len(df_plot)}")
    logging.info(pivot)

logging.info("Fetching variants with guide hits...")
variants_with_hits = gql(GET_VARIANTS_WITH_HITS, {"cds_group": CDS_GROUP})["variants"]
logging.info(f"✅ Found {len(variants_with_hits)} variants")


for p in [85,90,95,100]:
    plot_hits_consequences(variants_with_hits, guide_type="pre_mRNA", hits_threshold=p)

logging.info("Fetching variants with guide hits...")
variants_with_hits = gql(GET_VARIANTS_WITH_HITS, {"cds_group": CDS_GROUP})["variants"]
logging.info(f"✅ Found {len(variants_with_hits)} variants")

# Filter to only variants with mature_mRNA guides
variants_with_mature_mrna = [
    v for v in variants_with_hits 
    if any(vg["guide"]["type"] == "mature_mRNA" for vg in v.get("variants_guides", []))
]
logging.info(f"✅ Filtered to {len(variants_with_mature_mrna)} variants with mature_mRNA guides")

for p in [85,90,95,100]:
    plot_hits_consequences(variants_with_mature_mrna, guide_type="mature_mRNA", hits_threshold=p)



# ── Fetch good variants with bystander details ──
# Same filters as the off-target plot (GET_VARIANTS_WITH_HITS), just fetching bystander data

GET_GOOD_VARIANTS_BYSTANDERS = """
query GetGoodVariantsBystanders($cds_group: [String!]){
  variants(
    where: {
      class: {_eq: "SNV"}
      _or: [
        {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
        {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
      ]
      variants_features: {
        feature: {biotype: {name: {_eq: "protein_coding"}}}
      }
      gene: {
        genes_quantitative_scores: {
          quantitative_score: {name: {_eq: "SFARI Gene Score"}}
        }
      }
      variants_consequences: {
        consequence: {name: {_in: $cds_group}}
      }
    }
  ) {
    id
    variants_consequences {
      consequence { name }
    }
    variants_guides(where: {guide: {type: {_eq: "pre_mRNA"}}}) {
      bystanders {
        cadd_phred
      }
    }
  }
}
"""

good_variants = gql(GET_GOOD_VARIANTS_BYSTANDERS, {"cds_group": CDS_GROUP})["variants"]
logging.info(f"✅ Found {len(good_variants)} good variants (matching plot)")

# Check variants without guides
no_guide_count = sum(1 for v in variants_with_hits if not v.get("variants_guides", []))
logging.info(f"Variants without guides: {no_guide_count}")

# Check variants with NULL hits_85
null_hits_count = sum(1 for v in variants_with_hits 
                      if v.get("variants_guides", []) and 
                      v["variants_guides"][0]["guide"].get("hits_85") is None)
logging.info(f"Variants with NULL hits_85: {null_hits_count}")

# print 5 ids of variants without guides
logging.info("Example variant IDs without guides:")
for v in variants_with_hits:
    if not v.get("variants_guides", []):
        logging.info(f"  {v['id']}")


# ── Bystander burden plot function ──
from collections import Counter

bin_colors = {
    "0": colors['forestgreen'],
    "1": colors['yellowgreen'],
    "2": '#FFEA00',
    "3": colors['orange'],
    "4": colors['red2'],
    "5": colors['red2'],
    "6": colors['red2'],
    "7": colors['red2'],
    "8": colors['red2'],
    "9": colors['red2'],
    "10+": colors['red2'],
}

def plot_bystander_burden(variants: list[dict], cadd_cutoff: int, output_dir: str):
    """Plot bystander editing burden histogram for a given CADD phred cutoff."""
    bystander_counts = []
    for v in variants:
        variants_guides = v.get("variants_guides", [])
        if not variants_guides:
            continue
        bystanders = variants_guides[0].get("bystanders", [])
        n_harmful = sum(1 for b in bystanders if b.get("cadd_phred") is not None and b["cadd_phred"] >= cadd_cutoff)
        bystander_counts.append(n_harmful)

    bin_order = [f"{i}" for i in range(10)] + ["10+"]
    binned = [str(c) if c < 10 else "10+" for c in bystander_counts]

    count_per_bin = Counter(binned)
    values = [count_per_bin.get(b, 0) for b in bin_order]
    bar_colors = [bin_colors[b] for b in bin_order]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(bin_order, values, color=bar_colors, edgecolor='black', linewidth=0.8)

    for i, (bin_label, val) in enumerate(zip(bin_order, values)):
        if val > 0:
            ax.text(i, val + max(values) * 0.01, str(val), ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_xlabel(f"Number of bystanders with CADD phred ≥ {cadd_cutoff}", fontsize=12)
    ax.set_ylabel("Number of variants", fontsize=12)
    ax.set_title(f"Bystander Editing Burden for Editable Variants (CADD ≥ {cadd_cutoff})", fontsize=13)
    ax.set_ylim(0, max(values) * 1.08)
    plt.tight_layout()
    fig.savefig(f"{output_dir}/plots/bystander_burden_pre_mRNA_cadd{cadd_cutoff}.svg", format="svg", bbox_inches="tight")
    plt.show()

    logging.info(f"\nCADD ≥ {cadd_cutoff} — Total variants: {len(bystander_counts)}")
    logging.info(f"Breakdown: {dict(sorted(count_per_bin.items()))}")


for cutoff in [20, 25]:
    plot_bystander_burden(good_variants, cadd_cutoff=cutoff, output_dir=str(OUTPUT_DIR))

# Step 1: Get SFARI gene ENSGs (small, fast)
SFARI_GENES = """
query {
  genes(where: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}}) {
    ensg
  }
}
"""
sfari_ensgs = [g["ensg"] for g in gql(SFARI_GENES, {})["genes"]]
logging.info(f"SFARI genes: {len(sfari_ensgs)}")

# Step 2: Get CDS variants in those genes
GET_CDS_VARIANTS = """
query GetCdsVariants($cds_group: [String!], $ensgs: [String!]) {
  variants(
    where: {
      class: {_eq: "SNV"}
      gene_ensg: {_in: $ensgs}
      variants_features: {
        feature: {biotype: {name: {_eq: "protein_coding"}}}
      }
      variants_consequences: {
        consequence: {name: {_in: $cds_group}}
      }
    }
  ) {
    id
    variants_consequences {
      consequence { name }
    }
  }
}
"""

logging.info("Fetching CDS variants...")
cds_variants = gql(GET_CDS_VARIANTS, {"cds_group": CDS_GROUP, "ensgs": sfari_ensgs})["variants"]
logging.info(f"✅ Found {len(cds_variants)} variants")

QUERY_CDS_VS_OTHER_TOTAL = """
query CdsVsOtherTotal($cds: [String!]!) {
  total_variants: variants_aggregate {
    aggregate { count }
  }

  cds_group_variants: variants_consequences_aggregate(
    where: { consequence: { name: { _in: $cds } } }
  ) {
    aggregate {
      count(distinct: true, columns: variant_id)
    }
  }
}
"""
QUERY_CDS_VS_OTHER_PROTEIN_CODING_SFARI = """
query CdsVsOtherProteinCodingSfari($cds: [String!]!) {
  total_variants: variants_aggregate(
    where: {
    class : {_eq: "SNV"}
      variants_features: { feature: { biotype: { name: { _eq: "protein_coding" } } } }
      gene: {
        genes_quantitative_scores: {
          quantitative_score: { name: { _eq: "SFARI Gene Score" } }
        }
      }
    }
  ) {
    aggregate { count }
  }

  cds_group_variants: variants_consequences_aggregate(
    where: {
      consequence: { name: { _in: $cds } }
      variant: {
        class : {_eq: "SNV"}
        variants_features: { feature: { biotype: { name: { _eq: "protein_coding" } } } }
        gene: {
          genes_quantitative_scores: {
            quantitative_score: { name: { _eq: "SFARI Gene Score" } }
          }
        }
      }
    }
  ) {
    aggregate {
      count(distinct: true, columns: variant_id)
    }
  }
}
"""
CDS_GROUP = ["StopGained", "Missense", "Splice", "Frameshift", "Synonymous", "StartLost", "StopRetained"]

def run_and_print(title: str, query: str):
    res = gql(query, {"cds": CDS_GROUP})
    total = res["total_variants"]["aggregate"]["count"]
    cds = res["cds_group_variants"]["aggregate"]["count"]
    other = total - cds

    pct_cds = 100 * cds / total if total else 0
    pct_other = 100 * other / total if total else 0

    logging.info(f"\n{title}")
    logging.info(f"• CDS-group   : {cds:>10,}  ({pct_cds:5.2f}%)")
    logging.info(f"• Other       : {other:>10,}  ({pct_other:5.2f}%)")
    logging.info(f"• TOTAL       : {total:>10,}")

run_and_print("ALL variants — CDS-group vs Other", QUERY_CDS_VS_OTHER_TOTAL)
run_and_print("protein_coding + SFARI Gene Score — CDS-group vs Other", QUERY_CDS_VS_OTHER_PROTEIN_CODING_SFARI)


# ── Neighbor drill ("Improve" editing class), scored offline (SIFT, release 113) ──
# Populate neighbor_edits for non-G>A missense variants, then build Improve guides + bystanders
# for the variants whose best edit improves. Reproducible with the main pipeline — a fresh run
# re-creates them. Logic lives in the `neighbor_drill` package; the throwaway
# `neighbor_drill/populate_now.py` runs the same functions against an already-populated DB.
from neighbor_drill.neighbor_edits_pipeline import populate_neighbor_edits
from neighbor_drill.improve_guides_pipeline import populate_improve_guides
logging.info("=== Neighbor drill: neighbor_edits ===")
populate_neighbor_edits()
logging.info("=== Improve guides + bystanders ===")
populate_improve_guides()


