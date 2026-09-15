# Non-standard ALT notation is excluded, not repaired

**Decision.** Variants whose alternate allele is written in a non-standard form are dropped
during preprocessing rather than rewritten into a standard one.

VariCarta records some alternate alleles as `ref/N` or `N/ref` - the alternate nucleotide
written beside the reference - and a few with an `I:` prefix for insertions. Both contain
characters outside `[ATCG]`, so the standard-allele filter in `pipeline/phase1_preprocess.py`
removes them.

## Why not repair them

Repairing is defensible on its own terms: `REF=C, ALT=C/G` plainly means the alternate
allele is `G`, and a parser that understood the notation would keep the variant. Code that
did exactly this existed briefly.

It was written **after** the analysis was run. Verified from the repository's own history:

| Date | Event |
|---|---|
| 2026-05-09 04:43 | the published database is built |
| 2026-05-10 | `phase1_preprocess.py` first committed, with no repair |
| **2026-06-15** | commit `c7e82ef` introduces `_strip_slash_alt` |

The notebook code that was live when the database was built only *inspects* these rows
(`str(row["REF"]) in str(row["ALT"])`); it never rewrites them.

## What repairing would cost

Measured by rebuilding the database from raw inputs with the repair in place, then diffing
the variant sets. The published set is a **strict subset**: 84 variants gained, none lost.
Of the 93 raw rows at those positions, 90 have `/` in ALT and 2 an `I:` prefix.

| Quantity | Published | With repair |
|---|---|---|
| Unique variants | **329,279** | 329,363 |
| SNVs | 281,812 | 281,894 |
| Protein-coding | 190,709 | 190,791 |
| In SFARI genes | 34,292 | 34,374 |
| CDS + splice | **9,962** | 10,043 |
| Directly correctable | **2,421** | 2,439 |
| Correctable, no deleterious bystander | **891** | 899 |

Every delta is exactly the subset of those 84 surviving that filter, confirmed stage by
stage: 84 -> 82 SNVs -> 82 protein-coding -> 82 SFARI -> 81 CDS. There is no second cause.

**The three editing-class headlines are unaffected**: 894 read-through candidates, 681 with
no off-target, 404 with no deleterious bystander, 158 improvable missense variants, 201 and
1,417 individuals, and 4,644 distinct individuals overall all reproduce exactly either way.

## Consequence

The repository reproduces the published numbers. The 90-odd excluded variants are a known
limitation of the data loader, counted and logged at run time rather than silently dropped:

```
Excluded, ALT written as 'ref/N' or 'N/ref': ...
Excluded, ALT written with an 'I:' prefix:   ...
```

Anyone extending this analysis should decide deliberately whether to parse that notation.
Doing so is a change to the input set, not a bug fix, and it moves the funnel.
