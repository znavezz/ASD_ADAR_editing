# Decision 0002 — SIFT ≤0.05 is deleterious, >0.05 is tolerated

**Date:** 2026-08-02 · **Status:** ⚠️ **reversed 2026-08-26** · **Decided by:** the authors

## Reversal (2026-08-26)

**The published convention stands: deleterious is SIFT `<0.05`, tolerated is SIFT `≥0.05`.**
The authors chose to keep the convention the code actually ran, so the *Methods* sentence is
the one that changes — not the Results, not the figures, not the code.

Everything below is retained as the record of the reasoning; it is no longer the operative
decision. Concretely, this reversal means:

- **No number moves.** 158, 360, 201, 4,644 / 13.1% and the Figure 3B quadrants
  (158 / 1,009 / 608 / 148) are all published as-is and stay as-is.
- **No figure is regenerated.** The "Figure 3B must be re-plotted" consequence recorded below
  and in `manuscript_corrections.md` §6–§7 does not apply.
- **The oracle baseline is already correct** — it pins 158 / 360 / 201 / 4,644, which is now
  both the published *and* the adopted convention.
- **The Methods sentence is corrected instead**, in two places, from `≤0.05` to `<0.05` and
  from `≥0.05` to the tolerated side. See `manuscript_corrections.md`.

The substantive finding below is unaffected and still worth reporting: the headline is
insensitive to this choice (158 / 150 / 154 across the three defensible rules), and the real
error was that the Methods and the Results disagreed — which this reversal fixes from the
other side.

---

## The question

The Methods section and the analysis code disagreed about which side of 0.05 the boundary
falls on:

| | deleterious | tolerated |
|---|---|---|
| Methods, as drafted | SIFT **≤** 0.05 | SIFT > 0.05 |
| Results, as drafted | SIFT **<** 0.05 | SIFT ≥ 0.05 |
| The code, as run (`SiftScorer.benign_cutoff`, `neighbor_drill/core.py:91`) | SIFT **<** 0.05 | SIFT ≥ 0.05 |

This is not a cosmetic wording difference. **114 variants in the Improve universe have a SIFT
score of exactly 0.05** (750 across the whole dataset), so the boundary is load-bearing.

It matters twice, with *opposite* polarity, which is why the effect is not a simple increase:

1. **The variant filter** — "the variant's own residue is deleterious". Moving `<` to `≤`
   **loosens** it and lets more variants in.
2. **The edit scorer** — "the best edit is tolerated". Moving `≥` to `>` **tightens** it and
   pushes some edits out.

Applying only the first (the reading quoted during review) gives 179 and is **wrong** — it is
the threshold half-applied. Applied consistently to both, the headline *falls*.

## The decision

**Deleterious is SIFT ≤ 0.05; tolerated is SIFT > 0.05.** The Methods wording stands; the
Results wording and the code change to match it.

Rationale: `≤0.05` is the conventional SIFT cutoff as stated by the SIFT authors and used
throughout the literature, so the Methods sentence was the correct one. Aligning to the
convention costs eight variants and buys a threshold a reviewer will not query.

The authors will inform the co-authors of the resulting number changes.

## Downstream effect — every number, verified

Re-derived directly against the frozen analysis database, read-only. **Every published value
below reproduced exactly before the new convention was applied**, which is what makes the new
column trustworthy.

| Quantity | Published (`<` / `≥`) | Adopted (`≤` / `>`) |
|---|---|---|
| Improve headline — deleterious → tolerated | **158** ✅ | **150** |
| Improve tier 3 — best edit tolerated & better | **360** ✅ | **331** |
| Individuals carrying an improved variant | **201** ✅ | **188** |
| Distinct individuals, all three classes | **4,644** ✅ | **4,631** |
| …as a share of 35,581 | 13.05% → **13.1%** | 13.02% → **13.0%** |

Splice-excluded variants of the same queries (diagnostic, not published):

| Quantity | Published | Adopted |
|---|---|---|
| `ImproveAnalysisNoSplice` headline | **150** ✅ | **143** |
| `ImproveAnalysisNoSplice` tier 3 | **350** ✅ | **322** |
| `ImproveSubjectsNoSplice` | **193** ✅ | **181** |

**Unchanged** — these do not touch the threshold: 3,166 direct-repair individuals · 1,417
read-through individuals · 4,879 missense pool · 1,963 editable · and the entire funnel
329,279 → 281,812 → 190,709 → 34,292 → 9,962 → 2,650 / 2,421 / 1,092 / 894 / 681 / 404 / 891.

### How 158 becomes 150

The two directions nearly cancel, which is why the net move is small and why applying one
side alone is so misleading:

```
158  published
+ 21  variants whose OWN SIFT is exactly 0.05 — now deleterious, and their best edit improves
− 29  variants whose BEST EDIT scores exactly 0.05 — no longer counted as tolerated
= 150
```

### Figure 3B quadrants

The axes and the 0.05 gridlines are unchanged; only points sitting exactly on a line move.
The total is unchanged at 1,923 plottable variants.

| Quadrant | Published | Adopted |
|---|---|---|
| deleterious → tolerated | **158** ✅ | **150** |
| deleterious → deleterious | **1,009** ✅ | **1,071** |
| tolerated → tolerated | **608** ✅ | **556** |
| tolerated → deleterious | **148** ✅ | **146** |
| total | **1,923** ✅ | **1,923** |

If the figure prints its quadrant counts, it must be regenerated. If it only draws points and
gridlines, it is unaffected — no point moves, only its label.

## ⚠️ Two different 150s — do not conflate them

The new Improve headline is **150**. The published `ImproveAnalysisNoSplice` headline is
**also 150**. These are **different sets of variants** and the equality is a coincidence:

```
150  no-splice, published threshold
150  with-splice, adopted threshold
122  variants in both
 28  in each that the other does not contain
```

This is the same trap as the *CHD8* 40/43 in [0001](0001-consequence-precedence.md), where
two unrelated quantities happened to be 3 apart. If a future check finds "150" it must
establish *which* 150 before concluding anything.

## Why the boundary is genuinely ambiguous

Ensembl reports SIFT to **two decimal places**, so a stored 0.05 is any true value in
[0.045, 0.055). Both `SIFT=deleterious(0.05)` and `SIFT=tolerated(0.05)` appear in the raw VEP
output for this run — the annotator itself disagrees with itself at the boundary. Re-fetching
at full precision was investigated and is not possible: Ensembl does not serve more digits.

So neither convention is recoverably "correct" for the 114 boundary variants. The choice is
made on convention, not on evidence, and this file is the record of that.

## Consequences for the code

The boundary must not survive as a bare literal repeated across call sites:

- The comparison is defined **once**, in `figures/_helpers.R`, as `sift_tolerated()`. Every
  call site reads it, and the bin labels come from `sift_labels()` on the same switch, so a
  printed threshold can never disagree with the comparison that produced the counts. Before
  this, the comparison appeared in four places across three files.
- The convention is **selectable at run time** through `SIFT_TOLERATED_STRICT`, so both
  readings are reproducible from this repository rather than argued about:

```bash
Rscript figures/ver1/Figure_2/99_no_marginals.R                        # as published
SIFT_TOLERATED_STRICT=TRUE Rscript figures/ver1/Figure_2/99_no_marginals.R
```

- `queries/queries.txt` carries both readings as named operations
  (`ImproveAnalysis` / `ImproveAnalysis_SiftInclusive`), each annotated with its verified
  count, so the difference can be inspected in the data instead of inferred from prose.

## Verification

Read-only, against the frozen database. Writes are rejected by the server, not merely avoided:

```bash
docker exec -i -e PGOPTIONS='-c default_transaction_read_only=on' asd_adar_postgres \
  psql -U admin -d asd_adar_db -X -f - < <query>
```

The published-column values are additionally re-checked by:

```bash
python scripts/verify_counts.py verify
```

Note that the oracle baseline still records **158 / 360 / 201 / 4,644** — it pins the analysis
*as published*, and is deliberately **not** updated by this decision. Re-baselining is a
separate, explicit act, and would only happen if the analysis itself were re-run under the new
convention.

## Neither side is a mistake — and there is a third option

*Added 2026-08-16, and it bears on whether this decision should be taken at all.*

The question "is `<0.05` wrong?" has a definite answer: **no**. At a printed score of exactly
0.05, VEP's own output disagrees with itself:

| VEP's label at `SIFT=…(0.05)` | occurrences |
|---|---|
| `deleterious` | 353 |
| `deleterious_low_confidence` | 44 |
| `tolerated` | 317 |
| `tolerated_low_confidence` | 36 |

**397 deleterious and 353 tolerated, at the same printed number.** VEP is not inconsistent —
it decided from the full-precision value and then printed two decimals, and the printing threw
the distinction away. So `<0.05` and `≤0.05` are two conventions applied to a value that no
longer contains the answer, and neither recovers it.

What *was* an error is that the Methods and the Results said different things. That has to be
fixed whichever convention wins.

### The third option: use VEP's own label

The label is stored per variant (`qualitative_scores`) and per edit
(`neighbor_edit_scores.label`), so the analysis can ask VEP what it decided instead of
re-deciding from a rounded number:

| rule | Improve headline | what it rests on |
|---|---|---|
| `<0.05` deleterious, `≥0.05` tolerated — **published** | **158** | a cutoff we chose, on a value rounded to 2 dp |
| `≤0.05` deleterious, `>0.05` tolerated — **adopted here** | **150** | the same cutoff, other side |
| **VEP's own `deleterious` / `tolerated` label** | **154** | VEP's decision, made before rounding |

All three sit within eight variants of each other, which is itself worth reporting: the
headline is not sensitive to this choice in any way that changes the paper's claim.

The third is the only one that needs no convention from us, and it is the one to prefer if the
question is *which is most defensible* rather than *which is conventional*. Its one remaining
judgement is explicit and small: whether `*_low_confidence` counts, and the 154 above says it
does.

**Not adopted here, because it is the authors' call**, and because switching to labels changes
what the Methods sentence must say — it would describe a source's classification rather than a
threshold, which is a different claim about the method.

## Why the headline falls when the filter is loosened

*Added 2026-08-16.* The obvious objection: if `<0.05` becomes `≤0.05`, more variants qualify as
deleterious — so how does the count go *down*? Decomposed against the database, and the parts
reconcile exactly:

| | |
|---|---|
| published, `<0.05` deleterious and `≥0.05` tolerated | **158** ✅ reproduces |
| **lost** — already counted, but the winning edit scored *exactly* 0.05, so it is no longer tolerated | **−29** |
| **gained** — the variant's own SIFT is *exactly* 0.05, so it is now deleterious, and its best edit clears 0.05 outright | **+21** |
| adopted, `≤0.05` deleterious and `>0.05` tolerated | **150** ✅ |

158 − 29 + 21 = 150.

**Why losses outweigh gains**, which is the part worth understanding rather than tabulating:

* A **gain** needs a conjunction — the variant must sit exactly on the boundary *and* have an
  edit that clears it strictly.
* A **loss** needs only that the *winning* edit landed on the boundary. And the winner is
  chosen as the best-scoring option, so for a marginal variant it is frequently the lowest
  score that still counted as tolerated — which is exactly 0.05.

The boundary value is over-represented among winners *because* it was the cheapest way to
qualify. Tightening the tolerated side therefore bites harder than loosening the deleterious
side, and the two do not cancel.

SIFT is reported to two decimals, so "0.05" is really the bin `[0.045, 0.055)`. That is why the
boundary holds 750 variants at all, and why which side it falls on is a decision rather than a
rounding detail.

## What it also touches, beyond the counts

*Added 2026-08-16.* The table above lists counts. Two more consequences were checked, and one
of them needs work that regenerating a number will not cover.

**No published count depends on the rescue SIFT.** `rescue_sift` appears in exactly two
queries, `NonG2A_StopGained_Pass_OLD` (81) and `NonG2A_StopGained_Pass_Subjects_OLD` (129), and
both are superseded `_OLD` variants. The published `NonG2A_StopGained_Pass_Subjects` (1,417)
filters on consequence, biotype, SFARI and NMD only. So the Rescue funnel is unaffected even
though 142 variant-guide rows carry a rescue SIFT of exactly 0.05.

**Figure 3B does change, and it prints the threshold.** The panel highlights points with
`cadd_phred >= 25 & rescue_sift >= 0.05`, labels its bins `"SIFT≥0.05"` / `"SIFT<0.05"`, and
its caption reads *"n = … rows (variant × guide); … pass CADD ≥ 25 AND SIFT ≥ 0.05"* with the
count and percentage interpolated. Under the adopted convention:

| | published (`≥0.05`) | adopted (`>0.05`) |
|---|---|---|
| variant × guide rows passing CADD ≥ 25 and tolerated | **922 (9.4%)** | **786 (8.0%)** |

So the figure needs regenerating, not just the text: the highlighted points move, the caption's
count and percentage change, and the printed threshold in the bin labels and caption has to
change with them. That panel is `figures/ver1/Figure_2/00_panel_B_cadd_sift.rmd`.

**This panel is not in the manuscript.** It is retained as part of the ver-1 record; the
published Figure 3B is the missense-optimization scatter
(`figures/Figure_3/03_panel_B_missense_scatter.R`), whose SIFT comparison does route through
`sift_tolerated()` and therefore honours the switch above.
