"""
neighbor_drill — the "Improve" ADAR editing class.

Layers:
  * core.py                    — pure logic (enumerate A→G codon edits, score, pick best). Portable:
                                 no DB, no VEP; you inject the scores.
  * offline_vep.py             — offline VEP scoring + GRCh37 CDS→genomic mapping (needs a VEP cache
                                 + the release-87 GTF + genome FASTA; no DB).
  * neighbor_edits_pipeline.py — populate_neighbor_edits(): query variants → core → offline_vep → insert.
  * improve_guides_pipeline.py — populate_improve_guides(): guides + bystanders (offline SIFT) for the
                                 variants whose best edit improves.
  * populate_now.py            — throwaway backfill script for the already-populated DB (delete later).

phase3 imports the two *_pipeline functions so a fresh pipeline run reproduces everything.
"""

from .core import (
    GENETIC_CODE_3LETTER,
    translate,
    Scorer,
    SiftScorer,
    Score,
    EditOption,
    editable_a_positions,
    enumerate_neighbor_edits,
    codon_cds_start,
    build_edit_hgvs,
    set_score,
    evaluate_options,
)

__all__ = [
    "GENETIC_CODE_3LETTER", "translate", "Scorer", "SiftScorer", "Score", "EditOption",
    "editable_a_positions", "enumerate_neighbor_edits", "codon_cds_start", "build_edit_hgvs",
    "set_score", "evaluate_options",
]
