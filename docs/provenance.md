# Provenance of the published analysis

Exact versions, digests and file identities behind the numbers in `tests/golden/oracle/`.
Captured 2026-08-01 from the artifacts of the run itself, not from memory or documentation.

The frozen tree that produced these results is tagged **`paper-2026.0`**.

Most of this was recoverable only because VEP writes its own version block into its output.
Nothing in the pipeline recorded it deliberately — which is the gap the run manifest and the
generated Methods section are meant to close.

---

## Annotation — Ensembl VEP

Authoritative source: the `##` header of `Output/VEP/varicarta_vepped.txt`, written by the run.

| Item | Value |
|---|---|
| VEP | **v113.0** |
| Run at | 2026-05-07 20:03:48 |
| Cache | `homo_sapiens/113_GRCh37` |
| API | 113 |
| ensembl | 113.58650ec |
| ensembl-variation | 113.ebfce74 |
| ensembl-funcgen | 113.e30608c |
| ensembl-compara | 113.9cf749d |
| ensembl-io | 113.bee6816 |
| Assembly | **GRCh37.p13** |
| gencode | GENCODE 19 |
| genebuild | 2011-04 |

### Databases bundled in that cache

These are the versions actually used. The manuscript cites the method paper for each tool, as
is conventional; the column on the right notes where the *release* is also worth stating.

| Source | Version | In the manuscript |
|---|---|---|
| **ClinVar** | **202306** (June 2023) | cited (ref 30, Landrum 2018); release **not stated — worth adding** |
| SIFT | 5.2.2 | cited (ref 32, Ng & Henikoff 2003); version not stated (conventional) |
| PolyPhen | 2.2.2 | not used in the reported analyses |
| dbSNP | 156 | cited (ref 47); build not stated (conventional). **Genuinely used** — see below |
| 1000 Genomes | phase3 | — |
| COSMIC | 98 | not used |
| HGMD-PUBLIC | 20204 | not used |
| regbuild | 1.0 | not used |

**Only ClinVar is worth adding to the Methods.** ClinVar reclassifies variants substantially
between releases, and the manuscript makes counting claims about pathogenic / likely-pathogenic
variants (Figure 2B, and item 5 of the audit). Those counts are pinned to the June 2023
snapshot frozen into the Ensembl 113 cache, not to ClinVar at the time of writing. Tool
versions for SIFT and dbSNP are cited by method paper in the normal way and need no change.

Ensembl 113 itself **is** stated in the Methods and cited (ref 44, Dyer 2025) — correct as
written.

### Container

| Item | Value |
|---|---|
| Image | `ensemblorg/ensembl-vep` — **invoked untagged** (`phase1_preprocess.py:138`), i.e. `:latest` |
| Digest as cached locally | `sha256:f864f0a9ea30d77ea50f1a570fbf7ed7bd037bdd7da8637e0ba6228027639ecf` |
| Image created | 2024-11-25 |
| Local image id | `sha256:1d2b74402bf762db911b537c16df1307fd7fc70dcb5935090e54f582b6e65ef7` |

The image on disk predates the 2026-05-07 VEP run and `docker run` does not re-pull when an
image is present, so **the cached digest is the one that produced the published output** and
the link is sound. It is nonetheless recoverable only by inference. Pin it explicitly:

```
ensemblorg/ensembl-vep@sha256:f864f0a9ea30d77ea50f1a570fbf7ed7bd037bdd7da8637e0ba6228027639ecf
```

**The CADD version is properly cited, and now independently verified.** The manuscript states
"CADD v1.7" and cites the CADD v1.7 paper (ref 31, Schubach 2024). What was missing was the
*evidence trail*: CADD came from the Ensembl REST API (`?CADD=1`, `phase2_db.py:2400`), which
returns no version string, and nothing stored one — so the pipeline could not confirm which
build was served.

**Closed by direct comparison (2026-08-09).** Every stored CADD PHRED value was looked up in
the local `whole_genome_SNVs.tsv.gz` (CADD v1.7 GRCh37) and compared:

| | |
|---|---|
| stored `CADD_phred` values compared | **281,812** |
| found in the local v1.7 file | 281,812 (none missing) |
| identical to 3 decimal places | **281,812 — 100.0000%** |
| differing | **0** |

So the REST endpoint was serving CADD v1.7 GRCh37. The published claim is correct, and the
values are now reproducible from a checksummed local file instead of a live endpoint whose
build is unrecorded and will eventually change. Switching the pipeline to the local file
changes no number; it only makes the claim checkable.

The local file that was compared against is identified exactly:

| | |
|---|---|
| Path | `~aluguest/.vep/Plugins/CADD/v1.7_GRCh37/whole_genome_SNVs.tsv.gz` |
| Size | 85,228,947,819 bytes |
| MD5 (payload) | `fa3ee79df2f509ead49438f60209566d` ✅ matches upstream |
| MD5 (`.tbi`) | `6bc06b08fd9a521d05d970a06939e52d` ✅ matches upstream |
| Published at | `krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/MD5SUMs`, copied to `Resources/checksums/` |

Both digests were cross-checked against a second mirror (`kircherlab.bihealth.org`), which
agrees on these two files. It does **not** agree on the gnomAD indel tables in the same
directory — irrelevant here, since this analysis is SNV-only and never reads them.

> ⚠️ **A truncated decoy sits one directory up.**
> `~aluguest/.vep/Plugins/whole_genome_SNVs.tsv.gz` — same filename, 1.37 GB, no `.tbi`, no
> BGZF end-of-file marker. It stops mid-record at chr1:46,933,893, about 1.5% of the genome.
> Nothing in this repository reads it, and the VEP invocation loads AlphaMissense, Blosum62,
> LoFtool, MaveDB, NMD and pLI but **not** CADD — so no published number is affected.
> It matters for the rewrite: a tabix query against a truncated file returns *no rows* rather
> than failing, so pointing at it yields a complete-looking run whose scores are silently
> absent past that point, biased by genomic position. Pin the path *and* the checksum.

**Allele frequencies do not come from CADD or from gnomAD directly.** They arrive in VEP's
`Extra` column via `--af_gnomade`/`--af_gnomadg`/`--max_af` (`allele_freq_common.py:13`). On
the GRCh37 cache the gnomAD side is **exomes r2.1 only** — every `gnomadg_*` (genomes) column
is structurally NULL — and `MAX_AF` is a **VEP-computed** maximum across 1000G/ESP/gnomAD,
not a gnomAD-provided value.

### The rescue SIFT is the only score that needed the network

Two different SIFT numbers appear in this analysis, and only one of them came from REST:

* **the variant's own SIFT** — read from the offline VEP run's `Extra` column, like every
  other score. Never queried over the network.
* **the rescue SIFT** — the score of the read-through substitution (reference codon → `TGG`),
  which is a variant that does not exist in the input and therefore has no row in the VEP
  output. `_build_rescue_sift_query` (`phase3_guides.py:324`) builds a coding-HGVS string and
  queries `/vep/human/hgvs`.

**Why REST, and not the local cache.** Offline VEP refuses HGVS input outright — `ERROR:
Cannot use HGVS format in offline mode` — because resolving a transcript-relative coordinate
needs the database. So this was not a preference; the offline path could not answer the
question as asked. Doing it locally would mean converting the codon to genomic
coordinates first and submitting a VCF; that conversion is exactly what the comparison
below performs.

**Comparison (2026-08-09).** Each rescue codon was converted to genomic coordinates, verified
base-for-base against hg19, and re-scored with the same offline VEP 113 GRCh37 cache:

| | |
|---|---|
| rescue-SIFT variants compared | **4,864** of 4,897 |
| identical | **4,863 — 99.98%** |
| differing | **1** |
| not comparable | 33 (0.7%) — codon split across an exon boundary, so it has no contiguous genomic span |

The single disagreement is variant 297884 (`chr20:35,740,753`, *MROH8*, `ENST00000343811.4`):
REST 0.26, local 0.01 — which straddles the 0.05 cutoff. **It is outside every published
count**: *MROH8* carries no SFARI gene score, and `rescue_sift` is used in exactly two queries
(`queries.txt:281`, `:337`), both of which require SFARI membership.

So moving the rescue SIFT to the local cache changes no published number either. The 33 split
codons are the real implementation cost, and they need genuine exon structure — the `exons`
table in this database is empty.

---

## Alignment — BLAT

| Item | Value |
|---|---|
| Version | **Standalone BLAT v. 36x2** (matches the Methods) |
| Binary | resolved from `$BLAT_PATH` |
| Flags | `-out=psl -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0` (`phase3_guides.py:854-858`) |
| Identity thresholds | 85 / 90 / 95 / 100 %, computed by a Python port of UCSC `pslCalcMilliBad` |

`-minIdentity=0` means BLAT applies no identity filter itself; the 85/90/95/100 % cuts are
applied afterwards in pandas (`phase3_guides.py:918`). Identity is measured **over the aligned
block** — `pslCalcMilliBad` normalises by `matches + repMatches + misMatches`, not by guide
length — so a 20-nt perfectly-matching stretch of a 41-nt guide registers as a 100 %-identity
hit. This makes the off-target screen more permissive than a whole-guide identity cut would
be, and therefore conservative with respect to the "no off-target" counts.

The binary is not vendored: `$BLAT_PATH` points at a conda environment outside this
repository, so a fresh checkout must install BLAT 36x2 and set that variable. Pinning it as a
declared, checksummed dependency is an open reproducibility item.

---

## Input data

| Resource | File | Size |
|---|---|---|
| Variants | `varicarta_2024-12-12_11-46-40.vcf` | 463.3 MB |
| Reference genome | `all.fa` (hg19) | 3,157.6 MB |
| CDS reference | `Homo_sapiens.GRCh37.cds.all.fa` | 138.7 MB |
| MANE RNA | `MANE.GRCh38.v1.5.ensembl_rna.fna` | 78.1 MB |
| Gene annotation | `Homo_sapiens.GRCh37.87.gtf.gz` (release 87) | — |
| Gene scores | `SFARI-Gene_genes_07-08-2025release_07-17-2025export.csv` | — |
| Expression | GTEx v10 (`GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz`) | — |

Versions are encoded in filenames only; no checksums were recorded at the time of the run.
Checksums recovered since are collected in `Resources/checksums/`, alongside the URL each was
published at — a digest without a source proves a file is unchanged, not which release it is.

### Published checksums for every pinned input

*Added 2026-08-09.* A digest computed here proves a file has not changed since we looked at
it; only the **publisher's** digest proves which release it is. Everything upstream publishes
is now in `Resources/checksums/`, fetched verbatim.

| resource | upstream digest | verified |
|---|---|---|
| `Homo_sapiens.GRCh37.87.gtf.gz` | `Ensembl_GRCh37_release-87_gtf.CHECKSUMS` | ✅ `25146 41407` — matches |
| CADD v1.7 GRCh37 `whole_genome_SNVs.tsv.gz` + `.tbi` | `CADD_v1.7_GRCh37.MD5SUMs` | ✅ matched previously, two mirrors |
| VEP cache 113 GRCh37 | `Ensembl_release-113_vep_cache.CHECKSUMS` | ⚠️ published for the `.tar.gz`; the cache here is unpacked |
| `Homo_sapiens.GRCh37.cds.all.fa` | `Ensembl_GRCh37_release-87_fasta_cds.CHECKSUMS` | ⚠️ published for the `.gz`; the copy here is uncompressed |
| UCSC hg19 genome | `UCSC_hg19_bigZips.md5sum.txt` | ⚠️ published per archive; the copy here is a concatenated `all.fa` |
| `MANE.GRCh38.v1.5.ensembl_rna.fna` | — | ❌ NCBI publishes no digest for this release |
| SFARI gene export | — | ❌ per-query export, no published digest |
| VariCarta VCF | — | ❌ per-query export, no published digest |

Ensembl's `CHECKSUMS` files use BSD `sum`, not MD5: `sum <file>` reproduces both numbers.

**What ⚠️ means.** The digest is recorded and correct, but cannot be checked against the local
copy without re-downloading, because the local copy was unpacked or concatenated after
download. This is a real gap and not a hidden one: it is why the GTF, which is used as
distributed, is the only one of the four that could be verified outright.

### Assembly consistency of the guide searches

Both guide types are **extracted from GRCh37**. They are searched against different references:

| Guide | Cut from | Searched against | Consistent |
|---|---|---|---|
| `pre_mRNA` | `all.fa` (hg19 = GRCh37) | `all.fa` (GRCh37) | ✅ |
| `mature_mRNA` | `Homo_sapiens.GRCh37.cds.all.fa` (GRCh37) | `MANE.GRCh38.v1.5.ensembl_rna.fna` (**GRCh38**) | ❌ |

**The mismatch does not reach any published number.** Every query in `queries.txt` filters
`guide: {type: {_eq: "pre_mRNA"}}` — ten occurrences, and `mature_mRNA` appears in none of
them. The only figure code reading `mature_mRNA` hits is
`figures/ver1/Figure_1/03_panel_D_offtargets.rmd`, the **v1** panel, superseded by
`panels_D_E_offtargets_ver2.R`, which bins `dna_hits_85` (pre-mRNA).

So the affected quantity — `mature_mRNA` `hits_*` — is computed, stored, and never used in a
published count or figure. It is a defect in a dead branch, not in a result.

> An earlier version of this file said the mismatch was "carried into the published numbers".
> That was wrong: it was inferred from the mismatch existing, without checking which guide
> type the published queries actually filter on. Corrected 2026-08-03.

---

## Cohort sizes

| | records | distinct subjects |
|---|---|---|
| VariCarta source file | 525,687 | 35,581 |
| After preprocessing (staging table as loaded) | 525,516 | 35,492 |
| `subjects` table | — | 35,492 |

Preprocessing removed 171 records failing reference-allele validation, and 89 individuals with
them, before the staging load. The manuscript correctly reports the VariCarta figures.

---

## Coordinates

**1-based throughout, inclusive.** VCF `POS`, `coordinates.start`/`end`, VEP `cds_position`,
and `variants_guides.variant_idx` (the target sits at index 21 of a 41-nt guide) are all
1-based. 0-based indexing appears only inside pyfaidx slice expressions, always via an
explicit `-1` (`helpers.py:385`, `phase3_guides.py:182`). No BED, no pysam.

**Strand:** VEP reports alleles on the chromosome (+) strand. "Directly correctable by ADAR"
means G>A **on the coding strand**, i.e. `(G>A on + genes) OR (C>T on − genes)` in genomic
coordinates. Getting this wrong in one direction or the other is the cause of the 218-vs-138
discrepancy documented in `manuscript_number_audit.md`.

---

## Environment

- Python 3.10.20; exact package versions in `environment.lock.txt` (287 conda + 14 pip).
- Postgres 16, Hasura v2.40.0 (`docker-compose.yml`).
- R packages: not pinned; no `renv.lock` exists.

---

## Reproducing

```bash
python scripts/verify_counts.py verify    # 15 operations vs the committed baseline
python -m pytest                         # unit + golden tests, no network or database
```
