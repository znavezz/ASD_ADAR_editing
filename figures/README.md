# Figure code

R sources for the manuscript figures. The analysis lives in the Python pipeline and the
database; everything here reads results and draws.

## Running it

```bash
Rscript figures/run_all_figures.R
```

That rebuilds every figure end to end — panels first, then composites, then the supplement —
and needs no environment set: it derives the repository root from its own location. A run
prints what it rebuilt. It reports `[inherited]` if a panel came from a cached artifact rather
than being regenerated; a clean run reports none.

Set `REGEN=0` to reuse cached panel data where a panel supports it.

| Variable | What it is | Default |
|---|---|---|
| `ASD_PAPER_ROOT` | Repository root | derived from the script's location |
| `ASD_PAPER_OUT` | Where results are written | `$ASD_PAPER_ROOT/Output` |
| `ASD_ENV_FILE` | File holding `HASURA_ADMIN_SECRET` | `.env.local`, else `.env` |

Every path resolves through `_paths.R` and **fails loudly when unset** rather than writing
somewhere unexpected.

## Which directory makes which figure

Directory names state what they produce. The panel `.rds` filenames do **not** — they carry
the lettering of the layout they were first written for, and renaming them would orphan every
cached intermediate. Read this table, not the filenames.

| Directory | Produces |
|---|---|
| `Figures_1_2/` | **Figure 1** (mechanism illustration) and **Figure 2** (five data panels, A–E) |
| `Figure_3/` | **Figure 3** — nonsense rescue and missense optimization |
| `Figure_4/` | **Figure 4** — brain expression, treatability, CHD8 domains — and **Figure S1** |
| `ver1/` | Three panels the later scripts have no replacement for: the filtering funnel, the G>A consequence bar, and the brain-expression heatmap |
| `Lit_Search/` | The literature-audit pipeline behind Table 1 |

Figures 1 and 2 come from one script because they are one layout: the illustration sits above
the five data panels, and splitting them into two figures is a decision about presentation, not
about the analysis.

## Output

```
Output/Figures/
  Figure_1.{pdf,png,rds}   … Figure_4, Figure_S1
  _build/                  panels, per-panel PDFs, caches
```

`Figures/` holds one file per manuscript figure and nothing else. Everything needed to
reproduce is under `_build/` and `Output/Intermediate_tables/`.

## Reproducibility

Figures are **byte-identical between runs**. `ggrepel` places labels with random jitter, so
every call passes an explicit `seed`; without it the domain labels and their leader lines move
between runs while the data stays put. Nothing else in these scripts is stochastic.

The SIFT boundary is a convention rather than a measurement: scores are reported to two
decimal places, so a stored `0.05` is any true value in `[0.045, 0.055)`, and 750 variants sit
exactly on it.

**The published convention is the default, and needs no environment variable.** These scripts
treat `<0.05` as deleterious and `>=0.05` as tolerated, which is what the manuscript reports.

`SIFT_TOLERATED_STRICT=TRUE` renders the alternative `<=0.05`-deleterious reading. It is
provided so the difference can be inspected rather than argued about; **it is not the
published analysis** and no figure in the manuscript was produced with it set.
