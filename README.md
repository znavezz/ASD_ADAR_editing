# ADAR-mediated RNA base editing in Autism Spectrum Disorder

Analysis code for *"Exploring the Therapeutic Potential of ADAR-Mediated RNA Base Editors in
Autism Spectrum Disorder"*.

Variants from the [VariCarta](https://varicarta.msl.ubc.ca/) database are annotated, classified
by whether an ADAR A-to-I edit could correct or improve them, and assessed for off-target and
bystander editing risk. The repository holds the pipeline that produced every published number,
the queries that define each one, and the R code that draws every figure.

## Layout

| | |
|---|---|
| `pipeline/` | the analysis: three stages, plus `helpers.py` and the `neighbor_drill/` missense-optimization class |
| `scripts/` | the things you run — `run_pipeline.py`, `export_tables.py`, `verify_counts.py` |
| `queries/queries.txt` | **every published count, as a runnable query annotated with its verified value** |
| `figures/` | R sources for all figures — see `figures/README.md` |
| `db/schema.sql` | the database schema |
| `tests/` | the oracle baseline, a value-level subsample, and the query-parser tests |
| `Resources/manifest.tsv` | **every external input: release, size, checksum, source URL, and what attests it** |
| `Resources/checksums/` | the publishers' own checksum files, verbatim |
| `docs/provenance.md` | tool versions and data releases, recovered from the run's own artifacts |
| `docs/decisions/` | why the two conventions were chosen, with both readings costed |

Root holds only what has to live there: licence, environment and deployment config, and
pytest's own files. Everything else is a library (`pipeline/`), an entry point (`scripts/`),
data (`queries/`, `db/`, `Resources/`), or tests.

## Reproducing the published numbers

The analysis reads from a PostgreSQL database exposed through Hasura. With that stack running:

```bash
python scripts/verify_counts.py verify
```

This re-runs all 16 published counts and fails loudly if any disagrees with the committed
baseline. It never updates the baseline as a side effect; re-baselining is a separate,
explicit command.

The counts it checks are the paper's funnel and its three editing classes:

```
329,279 → 281,812 → 190,709 → 34,292 → 9,962      variant filtering
2,421 directly correctable; 891 bystander-clean; 1,923 (79%) with no off-target
1,092 → 894 nonsense → 681 no off-target → 404 no deleterious bystander
4,879 → 1,963 → 360 → 158 missense optimization
```

`queries/queries.txt` is the authority for what each number means. Every operation carries its verified
count as an annotation, so the file self-checks.

## Reproducing the figures

```bash
Rscript figures/run_all_figures.R
```

Rebuilds all five figures from the database. Output lands in `Output/Figures/`, one file per
manuscript figure. See `figures/README.md`.

## Environment

`environment.yml` declares the conda environment; `environment.lock.txt` records the exact
versions of all 287 conda and 14 pip packages used for the published run. Database credentials
are read from a local `.env` — see `.env.example`. No credential is committed.

## A note on two conventions

Two choices in this analysis are conventions rather than measurements, and both are documented
where they are made rather than buried:

- **The SIFT boundary.** Scores are reported to two decimal places, so a stored `0.05` is any
  true value in `[0.045, 0.055)`, and 750 variants sit exactly on it. The published analysis
  treats `<0.05` as deleterious. `queries/queries.txt` also carries the `≤0.05` reading, so both are
  reproducible from this repository and the difference can be inspected rather than argued.
- **Transcript choice.** Every published number is conditioned on VEP's `--pick` transcript.
  The alternative — reporting all transcripts — is a different analysis, not a correction.

## Authors

See `AUTHORS.md`.

## Citation

See the manuscript. `LICENSE` covers the code in this repository.
