"""
Offline VEP scoring + GRCh37 CDS→genomic mapping.

No DB and no project imports — reusable given: a VEP offline cache (release 113), the frozen
release-87 GRCh37 GTF (for CDS→genomic), and the genome FASTA. Offline VEP cannot take HGVS
input, so codon/base edits are expressed as **genomic** VCF records (mapped here), then scored
with `--sift b` from the cache.
"""

import os
import re
import gzip
import subprocess
from pathlib import Path

REVCOMP = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    return s.translate(REVCOMP)[::-1]


def to_dna(codon):
    """Codon string → DNA (uppercase, U→T)."""
    return codon.upper().replace("U", "T")


# ── GRCh37 CDS → genomic mapper (parse CDS features from the release-87 GTF) ──
def load_cds_intervals(gtf_path, wanted):
    """{transcript_id: {'strand': '+/-', 'cds': [(start, end), ...] genomic-ascending}} for `wanted`."""
    pid = re.compile(r'transcript_id "([^"]+)"')
    tx = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            c = line.split("\t")
            if c[2] != "CDS":
                continue
            m = pid.search(c[8])
            if not m or m.group(1) not in wanted:
                continue
            tx.setdefault(m.group(1), {"strand": c[6], "cds": []})["cds"].append((int(c[3]), int(c[4])))
    for t in tx.values():
        t["cds"].sort()
    return tx


def cds_to_genomic(t, pos):
    """1-based CDS coordinate → genomic position (handles split codons + strand)."""
    if t["strand"] == "+":
        rem = pos
        for s, e in t["cds"]:
            L = e - s + 1
            if rem <= L:
                return s + rem - 1
            rem -= L
    else:
        rem = pos
        for s, e in reversed(t["cds"]):
            L = e - s + 1
            if rem <= L:
                return e - rem + 1
            rem -= L
    return None


# ── offline VEP (SIFT) ──
def run_offline_vep(vcf_path, out_path, *, vep_cache, genome_fasta, work_dir, forks=8, assembly="GRCh37"):
    """Run VEP `--offline --cache --sift b` on a genomic VCF. `work_dir` must contain the VCF and
    is where the output is written (both are mounted into the container)."""
    vcf_path, out_path = Path(vcf_path), Path(out_path)
    vep_cache, genome_fasta, work_dir = Path(vep_cache), Path(genome_fasta), Path(work_dir)
    cmd = [
        "docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{vep_cache}:/data/vep",
        "-v", f"{work_dir}:/data/io",
        "-v", f"{genome_fasta.parent}:/data/fasta",
        "ensemblorg/ensembl-vep", "vep",
        "--offline", "--cache", "--dir_cache", "/data/vep",
        "--species", "homo_sapiens", "--assembly", assembly,
        "--fasta", f"/data/fasta/{genome_fasta.name}",
        "--fork", str(forks),
        "-i", f"/data/io/{vcf_path.name}", "-o", f"/data/io/{out_path.name}",
        "--sift", "b", "--transcript_version", "--force_overwrite", "--no_stats",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"VEP failed:\n{r.stderr[-2000:]}")


def parse_vep(path):
    """Parse a default-format VEP output file → {Uploaded_variation: [row_dicts]}."""
    res, header = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = line[1:].split("\t")
            continue
        if header is None:
            continue
        row = dict(zip(header, line.split("\t")))
        res.setdefault(row["Uploaded_variation"], []).append(row)
    return res


def row_for_transcript(rows, enst):
    """Pick the VEP output row for the variant's own transcript (falls back to any SIFT row)."""
    base = (enst or "").split(".")[0]
    for r in rows:
        if r.get("Feature", "").split(".")[0] == base:
            return r
    for r in rows:
        if "SIFT=" in r.get("Extra", ""):
            return r
    return rows[0] if rows else {}


def parse_sift(extra):
    """Extract (prediction, score) from a VEP `Extra` column, e.g. 'SIFT=tolerated(0.74)'."""
    for kv in (extra or "").split(";"):
        if kv.startswith("SIFT="):
            m = re.match(r"([a-z_]+)\(([\d.]+)\)", kv[5:])
            if m:
                return m.group(1), float(m.group(2))
    return None, None
