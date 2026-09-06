"""
Neighbor drill — the "Improve" editing class.

For a non-G>A **missense** variant, ADAR can edit the adenosines in the patient's
(alt) codon A→G. Each edit (or combination of edits) yields a different codon and a
different amino acid. We enumerate every option, let an injected, swappable *scorer*
(SIFT today, CADD/AlphaMissense tomorrow) judge each resulting amino acid, and pick
the best.

This module is intentionally **pure**: it has zero project imports (no DB, no env, no
VEP client). The pipeline passes plain inputs and an injected score-fetcher, so the
exact same logic runs inside phase3 today and as a standalone tool later.

Key conventions:
  * Codons are 3-letter strings; helpers accept DNA or RNA, any case (T is treated
    as U internally for translation).
  * Codon positions are 1-based (1, 2, 3).
  * An A→G edit of a *sense* codon can never produce a stop codon, so every option is
    always a scorable missense-or-synonymous change — never a truncation.
"""

from dataclasses import dataclass, field
from itertools import combinations
from typing import Protocol

__all__ = [
    "GENETIC_CODE_3LETTER",
    "translate",
    "Scorer",
    "SiftScorer",
    "Score",
    "EditOption",
    "editable_a_positions",
    "enumerate_neighbor_edits",
    "codon_cds_start",
    "build_edit_hgvs",
    "set_score",
    "evaluate_options",
]

# RNA codon (U) → 3-letter amino acid; stops are 'Ter'.
GENETIC_CODE_3LETTER = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Ter', 'UAG': 'Ter',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Ter', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}


def _rna(codon: str) -> str:
    """Uppercase + DNA→RNA (T→U) so a codon can be given as DNA or RNA, any case."""
    return codon.upper().replace("T", "U")


def translate(codon: str) -> str:
    """3-letter amino acid for a codon (accepts DNA or RNA, any case); 'Unknown' if invalid."""
    return GENETIC_CODE_3LETTER.get(_rna(codon), "Unknown")


# ── Scoring abstraction ──────────────────────────────────────────────────────
# A Scorer encapsulates everything score-specific: directionality (which way is
# "less damaging"), the benign threshold, and how to compute an improvement vs a
# baseline. Adding a new score later = one tiny Scorer subclass; nothing else changes.

class Scorer(Protocol):
    source: str
    def rank_value(self, v: float) -> float: ...                  # higher = LESS damaging (uniform sort key)
    def tolerated(self, v: float | None) -> bool: ...             # passes the benign threshold?
    def delta(self, value: float, baseline: float) -> float: ...  # signed; >0 = improvement over baseline


@dataclass
class SiftScorer(Scorer):
    """SIFT: score in [0,1]; higher = more tolerated; tolerated when >= 0.05."""
    source: str = "SIFT"
    benign_cutoff: float = 0.05

    def rank_value(self, v):  return v
    def tolerated(self, v):   return v is not None and v >= self.benign_cutoff
    def delta(self, v, base): return round(v - base, 4)

# e.g. a future, opposite-direction score:
#   @dataclass
#   class CaddScorer(Scorer):
#       source: str = "CADD"; deleterious_cutoff: float = 20.0
#       def rank_value(self, v):  return -v               # higher CADD = MORE damaging
#       def tolerated(self, v):   return v is not None and v < self.deleterious_cutoff
#       def delta(self, v, base): return round(base - v, 4)


@dataclass
class Score:
    source: str
    value: float | None = None
    label: str | None = None
    improvement: float | None = None    # signed; >0 = better than the pathogenic AA (per scorer.delta)


@dataclass
class EditOption:
    edited_positions: tuple[int, ...]   # 1-based codon positions edited A→G
    post_codon: str                     # RNA, e.g. "CGU"
    post_aa: str                        # 3-letter
    restores_reference: bool = False
    scores: dict[str, Score] = field(default_factory=dict)   # 'SIFT' -> Score, later 'CADD' -> Score, …
    is_best: bool = False
    improves: bool = False


# ── Enumeration ──────────────────────────────────────────────────────────────

def editable_a_positions(alt_codon: str) -> list[int]:
    """1-based positions of ALL adenosines in the alt (patient) codon.

    The variant's own position is editable too for non-G>A variants: A→G there yields
    G, which never restores the C/T reference, so it's a real Improve option, not a Fix.
    """
    c = _rna(alt_codon)
    return [i + 1 for i in range(3) if c[i] == "A"]


def enumerate_neighbor_edits(alt_codon: str, ref_aa: str) -> list[EditOption]:
    """Every A→G combination of the alt codon's adenosines.

    1–3 editable A's → 1, 3, or 7 options (2ᵏ−1). Each option carries its resulting
    codon + amino acid and whether it restores the reference amino acid (synonymous to
    wild-type). Scores are attached later.
    """
    c = _rna(alt_codon)
    pos = editable_a_positions(c)
    out: list[EditOption] = []
    for r in range(1, len(pos) + 1):
        for combo in combinations(pos, r):
            chars = list(c)
            for p in combo:
                chars[p - 1] = "G"
            post = "".join(chars)
            aa = translate(post)
            out.append(EditOption(combo, post, aa, restores_reference=(aa == ref_aa)))
    return out


# ── HGVS for scoring (mirrors phase3._build_rescue_sift_query) ────────────────

def codon_cds_start(cds_position, position_in_codon) -> int:
    """1-based CDS coordinate of the codon's first base (strand-agnostic: the stored
    codon is already coding-strand)."""
    return int(cds_position) - (int(position_in_codon) - 1)


def build_edit_hgvs(ref_codon: str, post_codon: str, cds_start: int, enst_ver: str) -> str:
    """Coding HGVS for the NET change reference→edited (single-sub or delins) so VEP
    scores the resulting amino acid relative to the reference. ref→post spans the
    variant position (ref→alt) plus the edited A→G position(s)."""
    r = _rna(ref_codon).replace("U", "T")
    p = _rna(post_codon).replace("U", "T")
    diffs = [i for i in range(3) if r[i] != p[i]]
    if len(diffs) == 1:
        i = diffs[0]
        return f"{enst_ver}:c.{cds_start + i}{r[i]}>{p[i]}"
    return f"{enst_ver}:c.{cds_start}_{cds_start + 2}delins{p}"


# ── Scoring + selection ──────────────────────────────────────────────────────

def set_score(option: EditOption, scorer: Scorer, value: float | None,
              label: str | None, baseline: float | None) -> None:
    """Attach a score to an option, computing its improvement vs the (per-source) baseline."""
    improvement = (scorer.delta(value, baseline)
                   if value is not None and baseline is not None else None)
    option.scores[scorer.source] = Score(scorer.source, value, label, improvement)


def evaluate_options(options: list[EditOption], primary: Scorer, variant_position: int) -> None:
    """Mutates options in place:

      (a) `improves` — a per-option fact: True for every restore-to-reference option and
          for every option whose primary score is tolerated AND beats the pathogenic AA.
      (b) `is_best` — the single winner. Restorers rank first; otherwise the best primary
          score. Restorers never short-circuit scoring the other options.

    Tie-break order: best score → fewest nucleotide changes → edits the variant's own
    position (the disease nucleotide).
    """
    def edits_variant(o: EditOption) -> bool:
        return variant_position in o.edited_positions   # True sorts ahead in max()

    for o in options:
        if o.restores_reference:
            o.improves = True                            # every WT-restoring option, however many
        else:
            s = o.scores.get(primary.source)
            if s is not None and s.value is not None:
                o.improves = primary.tolerated(s.value) and (s.improvement or 0) > 0

    restorers = [o for o in options if o.restores_reference]
    if restorers:                                        # Priority 1: restore-to-reference
        best = max(restorers, key=lambda o: (-len(o.edited_positions), edits_variant(o)))
    else:                                                # Priority 2: best primary score
        keyed = [o for o in options
                 if (s := o.scores.get(primary.source)) is not None and s.value is not None]
        if not keyed:
            return
        best = max(keyed, key=lambda o: (primary.rank_value(o.scores[primary.source].value),
                                         -len(o.edited_positions), edits_variant(o)))
    best.is_best = True
