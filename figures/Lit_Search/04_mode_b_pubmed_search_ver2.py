#!/usr/bin/env python3
"""
Ver-2 Mode B — broader PubMed search per top-20 gene.

Extends the ver-1 query set (8 templates) with 12 more to capture:
  - Missense-specific keywords (novel/pathogenic/rare variant)
  - Case-report contexts (WES/trio, DD/ID, epilepsy)
  - Loss-of-function language

Excludes PMIDs already covered by:
  - ver-1 triage corpus (per-variant .txt files)
  - ver-1 mode-B candidate list (ModeB_candidate_pmids.csv)

Output (ver-2 only):
  Results_Plots_Tables/Lit_Search/mode_b_ver2/ModeB_candidate_pmids_ver2.csv
  Intermediate_tables/Lit_Search/mode_b_ver2_raw_pmids_by_query.csv
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

TRIAGE_DIR    = VER1 / "Results_Plots_Tables/Lit_Search/Correctable/by_variant"
V1_CAND_CSV   = VER1 / "Results_Plots_Tables/Lit_Search/Correctable/ModeB_candidate_pmids.csv"
OUT_DIR       = VER2 / "Results_Plots_Tables/Lit_Search/mode_b_ver2"
RAW_DIR       = VER2 / "Intermediate_tables/Lit_Search"
OUT_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)

CANDIDATE_CSV = OUT_DIR / "ModeB_candidate_pmids_ver2.csv"
RAW_CSV       = RAW_DIR / "mode_b_ver2_raw_pmids_by_query.csv"

GENES = ["CHD8","ANK2","SCN2A","KMT2C","DSCAM","ASH1L","WDFY3","TSC2",
         "ZNF292","POGZ","ASXL3","CHD2","SHANK2","RELN","GRIN2B","MED13L",
         "ADNP","SHANK3","ARID1B","PTEN"]

# Ver-1 templates (kept for provenance) + 12 new — catches papers ver-1
# missed that describe rescue/missense/SIFT-Opt variants specifically.
QUERY_TEMPLATES = [
    # ver-1 originals
    "{g} nonsense",
    "{g} missense autism",
    "{g} splice variant",
    "{g} truncating variants",
    "{g} case report",
    "{g} variant patient",
    "{g} de novo mutation",
    "{g} functional analysis",
    # ver-2 additions
    "{g} pathogenic missense",
    "{g} novel variant",
    "{g} loss of function",
    "{g} premature stop codon",
    "{g} clinical exome sequencing",
    "{g} whole exome sequencing",
    "{g} developmental delay",
    "{g} intellectual disability",
    "{g} trio sequencing",
    "{g} pathogenic mutation",
    "{g} rare variant",
    "{g} likely pathogenic",
]
RETMAX = 50   # was 25 in ver-1

# --- 1. Existing PMID pool (must NOT re-audit these) ------------------------
PMID_RE = re.compile(r"PMID[: ]+([0-9]{5,9})")
existing_pmids = set()

for path in TRIAGE_DIR.glob("*.txt"):
    if ".failed." in path.name:
        continue
    existing_pmids.update(PMID_RE.findall(path.read_text(errors="replace")))
print(f"Ver-1 triage corpus PMIDs: {len(existing_pmids)}")

if V1_CAND_CSV.exists():
    with V1_CAND_CSV.open() as f:
        for row in csv.DictReader(f):
            if row.get("pmid"):
                existing_pmids.add(row["pmid"])
print(f"Existing PMID pool (triage + ver-1 mode B): {len(existing_pmids)}")

# --- 2. Esearch --------------------------------------------------------------
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

def esearch(query, retmax=RETMAX):
    params = {"db":"pubmed","term":query,"retmax":str(retmax),"retmode":"xml",
              "tool":"asd_mode_b_ver2","email":NCBI_EMAIL}
    url = f"{ESEARCH}?{urllib.parse.urlencode(params)}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            blob = r.read()
        root = ET.fromstring(blob)
        return [el.text for el in root.findall(".//IdList/Id") if el.text]
    except Exception as e:
        print(f"  query fail: {query} -> {e}")
        return []

# --- 3. Run all searches ----------------------------------------------------
raw_rows = []
new_pmids_per_gene = defaultdict(set)
n_searches = len(GENES) * len(QUERY_TEMPLATES)
print(f"Running {len(GENES)} genes × {len(QUERY_TEMPLATES)} queries = {n_searches} searches...")

for gi, gene in enumerate(GENES, 1):
    gene_new = set()
    for tmpl in QUERY_TEMPLATES:
        q = tmpl.format(g=gene)
        ids = esearch(q)
        for pmid in ids:
            raw_rows.append({"gene": gene, "query": q, "pmid": pmid,
                             "already_covered": "yes" if pmid in existing_pmids else "no"})
            if pmid not in existing_pmids:
                gene_new.add(pmid)
        time.sleep(0.34)
    new_pmids_per_gene[gene] = gene_new
    print(f"  [{gi:2}/20] {gene:<10} new PMIDs: {len(gene_new)}")

# --- 4. Write raw dump ------------------------------------------------------
with RAW_CSV.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","query","pmid","already_covered"])
    w.writeheader()
    w.writerows(raw_rows)

# --- 5. Distinct candidates + metadata fetch --------------------------------
all_new = sorted({p for s in new_pmids_per_gene.values() for p in s})
print(f"\nTotal distinct NEW candidate PMIDs across 20 genes: {len(all_new)}")

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 100
meta = {}
for i in range(0, len(all_new), BATCH):
    batch = all_new[i:i+BATCH]
    params = {"db":"pubmed","id":",".join(batch),"retmode":"xml","rettype":"abstract",
              "tool":"asd_mode_b_ver2","email":NCBI_EMAIL}
    url = f"{EFETCH}?{urllib.parse.urlencode(params)}"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            blob = r.read()
        root = ET.fromstring(blob)
        for art in root.iter("PubmedArticle"):
            pmid_el = art.find(".//PMID")
            if pmid_el is None:
                continue
            pmid = pmid_el.text.strip()
            ti = art.find(".//ArticleTitle")
            title = (ti.text or "")[:300] if ti is not None else ""
            authors = [a.find("LastName").text for a in art.findall(".//AuthorList/Author")
                       if a.find("LastName") is not None and a.find("LastName").text]
            year_el = art.find(".//Journal/JournalIssue/PubDate/Year")
            year = year_el.text if year_el is not None else ""
            jr_el = art.find(".//Journal/Title")
            jr = (jr_el.text or "") if jr_el is not None else ""
            meta[pmid] = {"title": title,
                          "first_author": authors[0] if authors else "",
                          "year": year, "journal": jr}
    except Exception as e:
        print(f"  meta fetch fail @{i}: {e}")
    time.sleep(0.34)

# --- 6. Write candidate table per (gene, PMID) ------------------------------
g_pmid_queries = defaultdict(set)
for r in raw_rows:
    if r["already_covered"] == "no":
        g_pmid_queries[(r["gene"], r["pmid"])].add(r["query"])

with CANDIDATE_CSV.open("w", newline="") as f:
    fields = ["gene","pmid","first_author","year","journal","title","queries_that_found_it"]
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    rows_count = 0
    for gene in GENES:
        for pmid in sorted(new_pmids_per_gene[gene]):
            m = meta.get(pmid, {})
            w.writerow({
                "gene": gene, "pmid": pmid,
                "first_author": m.get("first_author",""),
                "year": m.get("year",""), "journal": m.get("journal",""),
                "title": m.get("title",""),
                "queries_that_found_it": "; ".join(sorted(g_pmid_queries[(gene,pmid)])),
            })
            rows_count += 1

print(f"\nWrote {CANDIDATE_CSV}  rows={rows_count}")
print(f"Wrote {RAW_CSV}  rows={len(raw_rows)}")
print()
print(f"{'Gene':<10} {'new_PMIDs':>10}")
for g in GENES:
    print(f"{g:<10} {len(new_pmids_per_gene[g]):>10}")
