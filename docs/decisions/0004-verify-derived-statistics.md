# Decision 0004: statistics derived from a golden artifact are asserted, not just checksummed

**Date:** 2026-09-27 · **Status:** active · **Decided by:** the author

## Decision

`verify_counts.py verify` now recomputes the per-variant bystander statistics from
`Output/TableS1_Fixable_variants.csv` and fails if they move: 2,421 variants, 5.08 +/- 2.30
adenosines per variant, 1.35 +/- 1.47 predicted deleterious, and 891 variants with none.

## Why

The oracle covered counts returned by the GraphQL operations in `queries.txt`, and
`GOLDEN_ARTIFACTS` pinned the checksum of TableS1. Neither covered a number *computed from*
that table. The means and standard deviations above are quoted in the manuscript and in the
research proposal but were reachable only by rerunning an ad-hoc calculation, so a
transcription error would have survived both checks: the file hash still matches, and no
query returns the value.

## What it costs

The check reads a 2,421-row CSV on every `verify`, which is negligible, but it does mean
`verify` now fails if `Output/` is absent rather than silently skipping. That matches how
artifact verification already behaves.

Population standard deviation is used, matching how the values are quoted. Sample standard
deviation would give 2.2990 and 1.4697, which round identically at two decimals, so the
choice does not move a published number.

## Rejected

Adding an `avg`/`stddev` GraphQL operation to `queries.txt`, which the parser already
supports. It would have kept everything in one mechanism, but the statistic is over the
published table rather than over the database, and the table is the artifact a reader can
download. Checking the database would not catch the table and the text disagreeing.
