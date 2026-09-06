#!/usr/bin/env python3
"""
Ver-2 Mode B — generate per-(gene, PMID) Sonnet audit prompts.

Diffs from ver-1:
  - Reads the ver-2 candidate list (2,904 (gene, PMID) rows)
  - Reads the ver-2 631-variant top-20 treatable queue (RDS via rpy2? No —
    we use the CSV export written by 01_build_lit_queue_top20.R)
  - Caches PMC full-text at Intermediate_tables/Lit_Search/text_cache/<pmid>.txt
    so a re-run doesn't re-fetch (fetch is the slow step, ~30-45 min)
  - Writes prompts to Intermediate_tables/Lit_Search/mode_b_ver2/prompts/
"""
from __future__ import annotations
import csv
import re
import time
import os
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

# NCBI asks E-utilities callers to identify themselves. Supply a contact address
# via NCBI_EMAIL rather than committing one to a public repository.
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
if not NCBI_EMAIL:
    raise SystemExit("set NCBI_EMAIL to a contact address before querying PubMed")

VER1 = Path(os.environ["ASD_SOURCE"]) / "ver1"
VER2 = Path(os.environ["ASD_SOURCE"])

CAND_CSV    = VER2 / "Results_Plots_Tables/Lit_Search/mode_b_ver2/ModeB_candidate_pmids_ver2.csv"
QUEUE_CSV   = VER2 / "Intermediate_tables/Lit_Search/queue_top20_treatable.csv"
PROMPTS_DIR = VER2 / "Intermediate_tables/Lit_Search/mode_b_ver2/prompts"
INDEX_CSV   = VER2 / "Intermediate_tables/Lit_Search/mode_b_ver2/audit_index.csv"
TEXT_CACHE  = VER2 / "Intermediate_tables/Lit_Search/text_cache"

PROMPTS_DIR.mkdir(parents=True, exist_ok=True)
INDEX_CSV.parent.mkdir(parents=True, exist_ok=True)
TEXT_CACHE.mkdir(parents=True, exist_ok=True)

# Full-name AA → single-letter → three-letter
FULL_TO_ONE = {
    "Alanine":"A","Arginine":"R","Asparagine":"N","Aspartate":"D","Aspartic acid":"D",
    "Cysteine":"C","Glutamate":"E","Glutamic acid":"E","Glutamine":"Q","Glycine":"G",
    "Histidine":"H","Isoleucine":"I","Leucine":"L","Lysine":"K","Methionine":"M",
    "Phenylalanine":"F","Proline":"P","Serine":"S","Threonine":"T","Tryptophan":"W",
    "Tyrosine":"Y","Valine":"V","Stop":"*","Ter":"*"
}
ONE_TO_THREE = {"A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","E":"Glu","Q":"Gln",
                "G":"Gly","H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe",
                "P":"Pro","S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val","*":"Ter"}

GENES = {"CHD8","ANK2","SCN2A","KMT2C","DSCAM","ASH1L","WDFY3","TSC2","ZNF292",
         "POGZ","ASXL3","CHD2","SHANK2","RELN","GRIN2B","MED13L","ADNP","SHANK3",
         "ARID1B","PTEN"}

def to_one(name):
    if not name or name == "NA":
        return None
    if name in FULL_TO_ONE:
        return FULL_TO_ONE[name]
    return None

# --- Load the 631-variant queue and build variant identifiers per gene ------
gene_variants = defaultdict(list)
with QUEUE_CSV.open() as f:
    for r in csv.DictReader(f):
        gene = r["gene_symbol"]
        if gene not in GENES:
            continue
        ref1 = to_one(r.get("ref_aa_name",""))
        alt1 = to_one(r.get("alt_aa_name","")) if r.get("alt_aa_name","") not in ("","NA") else None
        pos  = r.get("protein_position","")
        chr_ = r.get("chr","")
        start= r.get("start","")
        tclass = r.get("treat_class","")
        is_nonsense = (tclass == "Nonsense Rescue")
        if is_nonsense and alt1 is None:
            alt1 = "*"
        # Stable variant_id used as prompt/output filename base
        if ref1 and alt1 and pos:
            variant_id = f"{gene}_p.{ref1}{pos}{alt1}_chr{chr_}_{start}"
        else:
            variant_id = f"{gene}_chr{chr_}_{start}"
        # HGVS strings
        ids = []
        if ref1 and pos:
            if is_nonsense:
                ids.append(f"p.{ref1}{pos}*")
                ids.append(f"{ref1}{pos}*")
                ids.append(f"p.{ref1}{pos}X")
                ids.append(f"{ref1}{pos}X")
                ref3 = ONE_TO_THREE.get(ref1,"")
                if ref3:
                    ids.append(f"p.{ref3}{pos}Ter")
                    ids.append(f"{ref3}{pos}Ter")
            elif alt1:
                ids.append(f"p.{ref1}{pos}{alt1}")
                ids.append(f"{ref1}{pos}{alt1}")
                ref3 = ONE_TO_THREE.get(ref1,"")
                alt3 = ONE_TO_THREE.get(alt1,"")
                if ref3 and alt3:
                    ids.append(f"p.{ref3}{pos}{alt3}")
                    ids.append(f"{ref3}{pos}{alt3}")
        ids.append(f"chr{chr_}:{start}")
        gene_variants[gene].append({
            "variant_id":     variant_id,
            "protein_change": r.get("hgvs_p","") or f"({tclass})",
            "consequence":    r.get("Consequence",""),
            "treat_class":    tclass,
            "ids":            ids,
        })

# --- Candidates -------------------------------------------------------------
candidates = []
with CAND_CSV.open() as f:
    for r in csv.DictReader(f):
        candidates.append(r)
distinct_pmids = sorted({c["pmid"] for c in candidates})
print(f"Candidates: {len(candidates)}  Distinct PMIDs: {len(distinct_pmids)}")

# --- Text fetch (with cache) -----------------------------------------------
ELINK  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def pmid_to_pmcid(pmid):
    params = {"dbfrom":"pubmed","db":"pmc","id":pmid,"retmode":"xml",
              "tool":"asd_mode_b_v2","email":NCBI_EMAIL}
    try:
        with urllib.request.urlopen(f"{ELINK}?{urllib.parse.urlencode(params)}", timeout=30) as r:
            blob = r.read()
        root = ET.fromstring(blob)
        for link in root.iter("Link"):
            pid = link.findtext("Id")
            if pid: return pid
    except Exception:
        pass
    return None

def fetch_pmc(pmcid):
    params = {"db":"pmc","id":pmcid,"retmode":"xml","rettype":"full",
              "tool":"asd_mode_b_v2","email":NCBI_EMAIL}
    try:
        with urllib.request.urlopen(f"{EFETCH}?{urllib.parse.urlencode(params)}", timeout=180) as r:
            blob = r.read()
        root = ET.fromstring(blob)
        return "\n".join(t for t in root.itertext())
    except Exception:
        return None

def fetch_abstract(pmid):
    params = {"db":"pubmed","id":pmid,"retmode":"xml","rettype":"abstract",
              "tool":"asd_mode_b_v2","email":NCBI_EMAIL}
    try:
        with urllib.request.urlopen(f"{EFETCH}?{urllib.parse.urlencode(params)}", timeout=60) as r:
            blob = r.read()
        root = ET.fromstring(blob)
        parts = []
        ti = root.find(".//ArticleTitle")
        if ti is not None and ti.text: parts.append("TITLE: " + ti.text)
        for at in root.findall(".//Abstract/AbstractText"):
            if at.text: parts.append(at.text)
        return "\n".join(parts)
    except Exception:
        return None

paper_texts = {}
n_pmc  = 0
n_abs  = 0
n_none = 0
print(f"Fetching text for {len(distinct_pmids)} PMIDs (cached to {TEXT_CACHE}) ...")
for i, pmid in enumerate(distinct_pmids, 1):
    cache_file = TEXT_CACHE / f"{pmid}.txt"
    source_file = TEXT_CACHE / f"{pmid}.source"
    if cache_file.exists() and cache_file.stat().st_size > 100:
        txt = cache_file.read_text(errors="replace")
        source = source_file.read_text().strip() if source_file.exists() else "cached"
        paper_texts[pmid] = (txt[:25000], source)
        if source.startswith("PMC"):
            n_pmc += 1
        elif source == "abstract":
            n_abs += 1
    else:
        pmc = pmid_to_pmcid(pmid)
        time.sleep(0.34)
        txt = None
        source = "unavailable"
        if pmc:
            txt = fetch_pmc(pmc)
            time.sleep(0.34)
            if txt and len(txt) > 200:
                source = f"PMC{pmc}"
                n_pmc += 1
        if not txt or len(txt) <= 200:
            abs_ = fetch_abstract(pmid)
            if abs_ and len(abs_) > 50:
                txt = abs_
                source = "abstract"
                n_abs += 1
            time.sleep(0.34)
        if txt:
            cache_file.write_text(txt)
            source_file.write_text(source)
            paper_texts[pmid] = (txt[:25000], source)
        else:
            n_none += 1
    if i % 100 == 0:
        print(f"  [{i}/{len(distinct_pmids)}]  pmc:{n_pmc} abs:{n_abs} none:{n_none}")

print(f"\nGot text for {len(paper_texts)} / {len(distinct_pmids)} PMIDs")
print(f"  pmc full-text: {n_pmc}   abstract-only: {n_abs}   no text: {n_none}")

# --- Prompt template --------------------------------------------------------
PROMPT_TEMPLATE = """TASK: Identify which of the listed variants (if any) are named in the paper text below. Be thorough — scan tables, prose, figure captions, and supplementary descriptions.

GENE: {gene}

VARIANTS of {gene} we want to check (one per line):
{variant_list}

PAPER:
  PMID:         {pmid}
  First author: {first_author}
  Year:         {year}
  Journal:      {journal}
  Source:       {source}

PAPER TEXT:
{text}

INSTRUCTIONS:
- For each variant listed above, decide if the paper names it in any form (HGVS short/long, rsID, chr:pos, residue#, prose).
- ALSO note any OTHER variants of {gene} mentioned in the paper that are NOT in our list (extract verbatim).
- If the paper isn't really about {gene}, mark UNRELATED.

OUTPUT FORMAT (exact):
  PAPER_RELEVANCE: <ABOUT_GENE | MENTIONS_GENE | UNRELATED>
  VARIANTS_NAMED_FROM_LIST: <variant_id1; variant_id2; ... or NONE>
  OTHER_GENE_VARIANTS_IN_PAPER: <verbatim snippets ; or NONE>
  JUSTIFICATION: <one line>
"""

# --- Write prompts ---------------------------------------------------------
written = 0
skipped = 0
with INDEX_CSV.open("w", newline="") as idxf:
    iw = csv.DictWriter(idxf, fieldnames=["gene","pmid","source","prompt_file","n_variants_listed"])
    iw.writeheader()
    for c in candidates:
        gene = c["gene"]
        pmid = c["pmid"]
        if gene not in gene_variants or not gene_variants[gene]:
            skipped += 1
            continue
        text, source = paper_texts.get(pmid, ("", "unavailable"))
        if not text:
            skipped += 1
            continue
        var_lines = []
        for vv in gene_variants[gene]:
            id_str = " | ".join(vv["ids"]) if vv["ids"] else "(no canonical IDs)"
            var_lines.append(f"  - {vv['variant_id']}  {vv['protein_change']}  [{id_str}]  [{vv['treat_class']}]")
        prompt = PROMPT_TEMPLATE.format(
            gene=gene, variant_list="\n".join(var_lines),
            pmid=pmid, first_author=c["first_author"],
            year=c["year"], journal=c["journal"], source=source,
            text=text,
        )
        out_path = PROMPTS_DIR / f"{gene}__{pmid}.prompt"
        out_path.write_text(prompt)
        iw.writerow({"gene": gene, "pmid": pmid, "source": source,
                     "prompt_file": out_path.name,
                     "n_variants_listed": len(gene_variants[gene])})
        written += 1

print(f"\nWrote {written} prompts to {PROMPTS_DIR}/")
print(f"Skipped {skipped} rows (missing gene variants or no text)")
print(f"Index: {INDEX_CSV}")
