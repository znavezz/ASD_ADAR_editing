"""Stage 1: read the VariCarta VCF, validate it against the genome, and annotate with VEP.

Turns the raw download into `varicarta_vepped.txt`, the single input every later stage
reads. Four things happen, in order:

  1. **Normalise and deduplicate.** 525,687 catalogued rows collapse to 329,411 unique
     (chr, pos, ref, alt) variants. Three rows where REF equals ALT are dropped: they
     describe no change.
  2. **Validate against hg19.** Each REF is checked against the reference base at that
     coordinate. 45 variants (0.01%) disagree and are removed rather than corrected -
     a REF that does not match the genome means the row's coordinates cannot be trusted.
  3. **Annotate.** VEP runs in a Docker container against the offline GRCh37 cache, with
     `--pick`, so exactly one transcript is chosen per variant. Every published number is
     conditioned on that choice.
  4. **Write** the annotated table for stage 2.

Coordinates are 1-based throughout, matching VCF. The only 0-based arithmetic is inside
pyfaidx slice expressions, always via an explicit `-1`.

Run through `scripts/run_pipeline.py`, which supplies the paths and the database target.
"""

# --- imports this stage needs -------------------------------------------------
# Previously inherited from the runner's scope via exec(). A stage that cannot say
# what it imports cannot be read, tested, or run on its own, so each now imports
# for itself; only genuine shared state still arrives from the runner.
import os
import re
import logging
import subprocess

df_varicarta = hp.load_vcf(VARICARTA_PATH)
# Strip whitespace from REF and ALT columns (some rows have trailing spaces)
df_varicarta["REF"] = df_varicarta["REF"].str.strip()
df_varicarta["ALT"] = df_varicarta["ALT"].str.strip()
_key_cols = ["#CHROM", "POS", "REF", "ALT"]

df_varicarta_unique = df_varicarta.drop_duplicates(subset=_key_cols)
logging.info(f"Total variants in Varicarta: {df_varicarta.shape[0]}")
logging.info(f"Unique variants by (#CHROM, POS, REF, ALT): {df_varicarta_unique.shape[0]} (dropped {df_varicarta.shape[0] - df_varicarta_unique.shape[0]} duplicates)")

# ── Normalize ALT field in df_varicarta (in place) ──
# Two non-standard ALT formats observed:
#   1. 'ref/N' or 'N/ref'  → keep only the alternate nucleotide N
#   2. 'I:NNN'             → strip the 'I:' prefix (insertion notation)
_slash_mask = df_varicarta_unique["ALT"].str.contains("/", na=False)
_insertion_mask = df_varicarta_unique["ALT"].str.startswith("I:", na=False)
logging.info(f"Variants with '/' in ALT (form 'ref/N' or 'N/ref'): {int(_slash_mask.sum())}")
logging.info(f"Variants with 'I:' prefix in ALT: {int(_insertion_mask.sum())}")

def _strip_slash_alt(row):
    parts = row["ALT"].split("/")
    if len(parts) == 2:
        if parts[0] == row["REF"]:
            return parts[1]
        if parts[1] == row["REF"]:
            return parts[0]
    return row["ALT"]

df_varicarta_unique.loc[_slash_mask, "ALT"] = df_varicarta_unique.loc[_slash_mask].apply(_strip_slash_alt, axis=1)
df_varicarta_unique.loc[_insertion_mask, "ALT"] = df_varicarta_unique.loc[_insertion_mask, "ALT"].str.slice(2)

# Duplicate counts: (#CHROM, POS, REF, ALT)
_dup_mask = df_varicarta_unique.duplicated(subset=_key_cols, keep=False)
logging.info(f"  ↳ of the '/' rows, duplicates by {_key_cols}: {int((_slash_mask & _dup_mask).sum())}")
logging.info(f"  ↳ of the 'I:' rows, duplicates by {_key_cols}: {int((_insertion_mask & _dup_mask).sum())}")

df_varicarta_unique = df_varicarta_unique.drop_duplicates(subset=_key_cols)
logging.info(f"Unique variants after ALT normalization: {df_varicarta_unique.shape[0]} (dropped {df_varicarta_unique.shape[0] - df_varicarta_unique.shape[0]} duplicates)")




pattern = re.compile(r'[^ATCG]')
df_varicarta_non_standard = df_varicarta_unique[df_varicarta_unique['ALT'].str.contains(pattern) | df_varicarta_unique['REF'].str.contains(pattern)]
logging.info(f"Non-standard Varicarta variants: {df_varicarta_non_standard.shape[0]}")

df_varicarta_standard = df_varicarta_unique[~(df_varicarta_unique['ALT'].str.contains(pattern) | df_varicarta_unique['REF'].str.contains(pattern))]
logging.info(f"Standard Varicarta variants: {df_varicarta_standard.shape[0]}")

# Remove variants where REF == ALT (not actual variants)
invalid_variants = df_varicarta_standard[
    df_varicarta_standard["REF"] == df_varicarta_standard["ALT"]
]
logging.info(f"⚠️  Found {len(invalid_variants)} invalid variants where REF == ALT:")
logging.info(invalid_variants[["#CHROM", "POS", "REF", "ALT", "ID"]])

logging.info(f"Total Varicarta variants: {df_varicarta.shape[0]}")
logging.info(f"Unique Varicarta variants: {df_varicarta_standard.shape[0]}")

# Filter them out
df_varicarta_standard_clean = df_varicarta_standard[
    df_varicarta_standard["REF"] != df_varicarta_standard["ALT"]
].copy()
logging.info(f"Num rows before removing invalid variants: {df_varicarta_standard.shape[0]}")
logging.info(f"Num rows after removing invalid variants: {df_varicarta_standard_clean.shape[0]}")
logging.info(f"Num dropped invalid variants: {df_varicarta_standard.shape[0] - df_varicarta_standard_clean.shape[0]}")





forks = int(os.environ["MAX_FORKS"])
logging.info(f"Using {forks} forks for multiprocessing.")   
df_varicarta_standard_clean_with_ref = df_varicarta_standard_clean.copy()
hp.add_genome_ref_column(
    df_varicarta_standard_clean_with_ref,
    fasta_path=GENOME_REFERENCE_HG19,
    genome_version="hg19",
    workers=forks
)
# Calculate percentage of mismatched variants
num_mismatches = len(df_varicarta_standard_clean_with_ref[df_varicarta_standard_clean_with_ref["REF"] != df_varicarta_standard_clean_with_ref["hg19"]])
total_variants = len(df_varicarta_standard_clean_with_ref)
mismatch_percentage = (num_mismatches / total_variants) * 100

logging.info(f"Mismatched variants: {num_mismatches}/{total_variants} ({mismatch_percentage:.2f}%)")


df_varicarta_standard_positive = df_varicarta_standard_clean_with_ref.copy()
df_varicarta_standard_positive.drop(df_varicarta_standard_positive[df_varicarta_standard_positive["REF"] != df_varicarta_standard_positive["hg19"]].index, inplace=True)
logging.info(f"Num rows before filtering: {df_varicarta_standard_clean_with_ref.shape[0]}")
logging.info(f"Num rows after filtering: {df_varicarta_standard_positive.shape[0]}")
logging.info(f"Num dropped rows: {df_varicarta_standard_clean_with_ref.shape[0] - df_varicarta_standard_positive.shape[0]}")


df_varicarta_standard_positive["Original_ID"] = df_varicarta_standard_positive["ID"]
df_varicarta_standard_positive["ID"] = df_varicarta_standard_positive["#CHROM"] + ":" + df_varicarta_standard_positive["POS"].astype(str) + ":" + df_varicarta_standard_positive["REF"] + ":" + df_varicarta_standard_positive["ALT"]
logging.info(f"Num unique variants: {df_varicarta_standard_positive['ID'].unique().shape[0]}")




path_varicarta_standard_positive = OUTPUT_DIR / "varicarta_standard_positive_clean_variants.vcf"
# Save to TSV without hg19 and Original_ID columns to ensure proper VCF format for VEP input
df_varicarta_standard_positive.drop(columns=["hg19", "Original_ID"]).to_csv(path_varicarta_standard_positive, sep='\t', index=False)



# ── VEP Annotation (replaces the %%bash cell) ──────────────────────────────
vep_output_file_name = "varicarta_vepped"

vep_output_dir = OUTPUT_DIR / "VEP"
vep_output_dir.mkdir(parents=True, exist_ok=True)

out_txt      = vep_output_dir / f"{vep_output_file_name}.txt"
stats_html   = vep_output_dir / f"{vep_output_file_name}.stats.html"
warnings_txt = vep_output_dir / f"{vep_output_file_name}.warnings.txt"


fasta_host = GENOME_REFERENCE_HG19

cache_host = VEP_CACHE


# ── Build the docker command ──
input_dir  = path_varicarta_standard_positive.resolve().parent
input_file = path_varicarta_standard_positive.name

uid_gid = f"{os.getuid()}:{os.getgid()}"

cmd = [
    "docker", "run", "--rm",
    "-u", uid_gid,
    "-v", f"{cache_host.resolve()}:/data/vep",
    "-v", f"{input_dir}:/data/input",
    "-v", f"{vep_output_dir}:/data/output",
    "-v", f"{fasta_host.resolve().parent}:/data/fasta",
    "ensemblorg/ensembl-vep", "vep",
        "--offline",
        "--cache", "--dir_cache", "/data/vep",
        "--fork", str(forks),
        "--buffer_size", "10000",
        "--species", "homo_sapiens", "--assembly", "GRCh37",
        "--fasta", f"/data/fasta/{fasta_host.name}",
        "-i", f"/data/input/{input_file}",
        "-o", f"/data/output/{out_txt.name}",
        "--stats_file", f"/data/output/{stats_html.name}",
        "--warning_file", f"/data/output/{warnings_txt.name}",
        "--pick",
        "--transcript_version",
        "--everything",
        "--sift", "b", "--polyphen", "b", "--protein",
        "--plugin", "AlphaMissense,file=/data/vep/Plugins/AlphaMissense_hg19.tsv.gz,cols=all",
        "--plugin", "Blosum62",
        "--plugin", "LoFtool,/data/vep/Plugins/LoFtool_scores.txt",
        "--plugin", "MaveDB,file=/data/vep/Plugins/MaveDB_variants.tsv.gz",
        "--plugin", "NMD",
        "--plugin", "pLI,/data/vep/Plugins/pLI_values.txt",
        "--custom", "file=/data/vep/homo_sapiens/113_GRCh37/custom/hg19.100way.phastCons.bw,short_name=PhastCons100,format=bigwig",
]

logging.info("Running VEP …")
result = subprocess.run(cmd, capture_output=True, text=True)

# Stream stderr (VEP progress) to log
if result.stderr:
    for line in result.stderr.splitlines():
        logging.info(f"[VEP] {line}")

if result.returncode != 0:
    raise RuntimeError(f"VEP failed with exit code {result.returncode}")

logging.info(f"✓ VEP finished: {out_txt.name}")
logging.info(f"✓ Stats file  : {stats_html.name}")


df_varicarta_annotated = hp.load_vcf(vep_output_dir / f"{vep_output_file_name}.txt", "Uploaded_variation")


logging.info(f"Number of annotated Variants: {df_varicarta_annotated.shape[0]}")
logging.info(f"PREPROCESS COMPLETED SUCCESSFULLY")
