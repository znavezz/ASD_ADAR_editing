# Decision 0001 — StopGained outranks Splice for dual-annotated variants

**Date:** 2026-08-02 · **Status:** accepted · **Decided by:** the authors

## The question

VEP can assign a variant more than one consequence. A variant at a splice site that also
creates a premature stop is annotated `StopGained,Splice`. Any figure or table showing
"variants by consequence" must place it in exactly one bucket, so a precedence rule is needed.

## The decision

**StopGained outranks Splice.** A variant annotated both is counted as **nonsense**.

Rationale: the therapeutic claim attached to these variants is stop-codon read-through, which
depends on there being a premature stop to recode. Classifying such a variant as "splice"
hides the property the analysis acts on.

## This ratifies existing behaviour — no published number changes

The analysis **already applies this rule, consistently, everywhere**. It is implemented once,
in the R prep, and every figure inherits it
(`figures/ver1/00_prep_data.rmd`):

```r
df_cds$New_Consequence <- dplyr::case_when(
  grepl("StopGained", df_cds$Consequence) ~ "Nonsense",   # tested FIRST
  grepl("Splice",     df_cds$Consequence) ~ "Splice",
  ...
)
```

Verified across the whole snapshot: all **36** `StopGained,Splice` variants carry
`New_Consequence == "Nonsense"`. There are no exceptions.

| Location | Value | Follows this rule? |
|---|---|---|
| Figure 2C — the 2,650 split | 1,994 / 333 / 228 / **94** / 1 | ✅ reproduces exactly |
| *CHD8* — "35 G>A, an additional 40 nonsense" | 35 + 40 + 81 = 156 | ✅ see below |
| Figure 4B / S1 gene bars | per-gene splits | ✅ same `New_Consequence` |

**The CHD8 numbers are correct.** *CHD8* has **43** variants bucketed as nonsense, of which
**3 are also G>A** and are counted under direct repair, since Figure 4B's precedence is G>A
first, then nonsense:

```
 35  G>A (Direct Repair)
 40  Nonsense Rescue     = 43 nonsense − 3 already counted as G>A
 81  Untreatable
156  total
```

The word "additional" in the sentence is doing real work: the 40 is the remainder after the
35, not a subset of the 43. Coincidentally *CHD8* also has exactly 3 `StopGained,Splice`
variants, but they are **different variants** — a numerical coincidence that made this look
like an inconsistency during review. It is not one.

## Two constants in the Python were out of step — now fixed

They turned out to be different situations, not two copies of one problem:

| File | Symbol | Was | Now |
|---|---|---|---|
| `pipeline/phase2_db.py:924` | `variants_rank` | `Splice` first, **defined but never read** | **deleted** |
| `scripts/export_tables.py` | `CONSEQUENCE_RANK` | `Splice` first | **StopGained first** |

`variants_rank` was dead: the only other mention in the repository was a comment in
`scripts/export_tables.py` saying it mirrored it. Removed rather than corrected.

`CONSEQUENCE_RANK` is live, but collapses **bystander** consequences, not variant ones — and
**no bystander in the dataset carries both Splice and StopGained** (0 of 359,291). So
reordering it is provably inert on this data while making the code state the rule the analysis
actually applies.

Both changes verified: the oracle passes unchanged.

In the rewrite this ordering exists **once**, as data rather than as a literal in two files.
One rule in two copies, free to drift, is the failure mode behind the strand convention (eight
copies) and the four `gql` variants (three different retry policies).

## How to justify differences later

A difference attributable to this rule must satisfy **all** of:

1. It is a *bucketed* count — variants shown in exactly one consequence category. Membership
   counts are unaffected, because `variants_consequences` matches if **any** consequence is
   StopGained. Verified: the 1,092 read-through pool contains 35 dual-annotated variants under
   either rule.
2. It moves nonsense up and splice down.
3. Its magnitude is within the dual-annotated bounds: **1** inside the 2,650, **3** in *CHD8*,
   **8** across the top-20 genes, **36** across the whole 9,962 CDS set.

Anything outside those bounds is a different problem.

**Unchanged by this decision:** 329,279 · 281,812 · 190,709 · 34,292 · 9,962 · 2,650 · 2,421 ·
1,092 · 894 · 681 · 404 · 4,879 · 1,963 · 360 · 158 · 891 · 1,923 — and the *CHD8* 35 / 40 / 81.

## Verification

```bash
python scripts/verify_counts.py verify
```
