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

This re-runs every annotated count and fails loudly if any disagrees with the committed
baseline: 18 values belonging to the numbers the manuscript reports, and 24 more that support
them, across 15 operations. It never updates the baseline as a side effect; re-baselining is a separate,
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

## What a full run needs, and how long it takes

The published numbers can be re-checked in seconds against an existing database
(`scripts/verify_counts.py verify`). Rebuilding that database from the raw inputs is a
different proposition, and these are the real costs.

**External data** (none of it is in this repository; `Resources/manifest.tsv` records each
one's release, checksum and source):

| | size |
|---|---|
| Ensembl VEP cache, `113_GRCh37` | **43 GB** |
| UCSC hg19 reference FASTA | 3.0 GB |
| GTEx v10 (TPM matrix + annotations) | 2.2 GB |
| VariCarta VCF | 442 MB |
| Ensembl GRCh37 release-87 CDS FASTA | 133 MB |
| MANE v1.5 transcripts | 75 MB |
| Ensembl GRCh37.87 GTF | 41 MB |
| SFARI gene scores | 136 KB |

**Software**: Docker (for PostgreSQL, Hasura and VEP), BLAT on `$PATH`, and the conda
environment in `environment.yml`. Network access is required: CADD and the read-through
SIFT scores are fetched from the Ensembl REST API.

**Disk for outputs**: about 3.5 GB, split between `Output/` (~900 MB) and the PostgreSQL
data directory (~2.4 GB).

**Measured wall-clock**, on 64 worker processes:

| stage | time |
|---|---|
| Phase 1, preprocess + VEP annotation | **~5 min** |
| Phase 2, database insertion | **~7.5 h** |
| Phase 3, guides, bystanders and BLAT off-targets | not recorded; dominated by BLAT |
| All annotated counts (`verify_counts.py`) | seconds |
| All five figures (`figures/run_all_figures.R`) | ~10 min |

Phase 2 is the long pole because every insert is a GraphQL mutation and CADD arrives from
the network in batches. Budget a day for a full rebuild, and note that it is restartable:
`--start_from db` and `--start_from guides` skip completed stages.

## Choosing which database a run writes to

`scripts/run_pipeline.py` and `reset_db.sh` both write, and both take `--env-file`
(default `.env`, or `$ASD_ENV_FILE`). The env file names the target database, so a
verification run can be pointed at a scratch stack without editing any source:

Both default to `.env.test`, the scratch stack, because both destroy or overwrite
data. Writing to the database named in `.env` is refused unless asked for by name:

```bash
./reset_db.sh --before 1            # resets the scratch stack
python scripts/run_pipeline.py      # rebuilds it

ASD_ALLOW_PUBLISHED_RESET=yes ./reset_db.sh --before 1 --env-file .env   # the real one
```

`--before 1` is a full teardown: it truncates the guide tables, deletes the PostgreSQL
data directory and the staging CSVs, and removes the VEP annotation. `--before 2` keeps
the VEP output; `--before 3` removes only guides and bystanders.

## Environment

`environment.yml` declares the conda environment; `environment.lock.txt` records the exact
versions of all 287 conda and 14 pip packages used for the published run. Database credentials
are read from a local `.env` — see `.env.example`. No credential is committed.

## A note on transcript choice

Every published number is conditioned on VEP's `--pick` transcript: one transcript per variant,
chosen by VEP. Reporting all transcripts is a different analysis, not a correction. This is a
convention rather than a measurement, so it is stated here rather than left implicit.

## Authors

See `AUTHORS.md`.

## Citation

See the manuscript. `LICENSE` covers the code in this repository.
