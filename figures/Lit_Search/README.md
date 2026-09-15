# Literature search

Table 1 of the manuscript lists 27 ADAR-correctable variants for which a publication
describes the variant in an affected individual. This directory holds the code that found
them and the evidence that substantiates each one.

## Where Table 1 comes from

```
evidence/Top20_named_papers_per_variant.csv   filtered to in_correctable == "yes"
  -> 27 variants across 10 genes
```

That file has one row per variant with `n_papers` and `pmids`. Its other 29 rows are
variants named in publications but outside the correctable set, so they are not in Table 1.

`evidence/Top20_named_papers_COMBINED.csv` is the long form: one row per
(variant, publication) pair, carrying `pmid`, `evidence_source`, `verification_method`,
first author, year, journal and title.

Verified 2026-09-15: the 27 genomic positions in that file match the 27 rows of the
manuscript's Table 1 exactly, in both directions.

## What is deterministic and what is not

The search universe is the 631 ADAR-correctable variants in the 20 most frequently affected
SFARI genes.

**Retrieval was LLM-assisted and is not bit-reproducible.** A Claude Sonnet 4.6 agent was run
per variant (May-June 2026) over structured inputs - coordinates, dbSNP identifiers, ClinVar
classification, in silico scores - and returned candidate publications. Re-running it would
return a different set: there is no seed, and the model version is not archivable.

**The variant-to-publication assignment is deterministic.** Retrieved PMIDs were checked
against NCBI records, then abstracts and, where available, full texts were searched for
variant-specific identifiers: HGVS nomenclature, rsIDs and genomic coordinates.
`03_extract_variants_from_mode_b.R` builds those patterns per variant and matches with
`grepl`, recording the pattern that fired. A citation was retained only where the variant is
explicitly described and linked to a human disease context.

So a reviewer can check every pair in Table 1 without running any model: the
`verification_method` column records which route each row came through (`regex`,
`PMC full-text`, `mode_a_sonnet`, `mode_b_sonnet`, `LLM (abstract)`), and the deterministic
matcher is the gate that decides membership.

## The scripts

| | |
|---|---|
| `01_build_lit_queue_top20.R` | build the queue: correctable variants in the top-20 SFARI-1 genes |
| `02_compare_ver1_coverage.R` | diff against the earlier triage table to isolate what needs new searching |
| `03_extract_variants_from_mode_b.R` | **the deterministic matcher** - builds HGVS/coordinate/rsID patterns and records which matched |
| `04_mode_b_pubmed_search_ver2.py` | PubMed E-utilities search. Set `NCBI_EMAIL` to a contact address first |
| `05_mode_b_audit_generate_ver2.py` | generate the per-(gene, PMID) extraction prompts |
| `06_mode_b_audit_runner_ver2.sh` | run them, concurrency 5, resumable |

They are not part of `run_all_figures.R`: the retrieval step is not reproducible, so it is not
wired into a command that claims to rebuild the paper. The evidence they produced is committed
here instead.
