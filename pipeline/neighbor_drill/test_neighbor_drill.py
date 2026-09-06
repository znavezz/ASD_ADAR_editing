"""
Unit tests for the pure neighbor_drill module (no DB / no network).

Standard-library `unittest` — no third-party dependency needed. Also collected by pytest
(see pytest.ini), so `pytest` from the repository root runs these alongside tests/.

Run **from the repository root** — this file imports the `neighbor_drill` package, so the
package's parent must be on sys.path:

    python -m unittest neighbor_drill.test_neighbor_drill -v
    python -m unittest discover
    python -m pytest neighbor_drill/test_neighbor_drill.py

Running from inside neighbor_drill/ fails: that puts the package's own directory on
sys.path rather than its parent, so `from neighbor_drill import ...` raises
ModuleNotFoundError. Requires Python >= 3.10 (core.py uses PEP 604 `X | None` in a class
body, with no `from __future__ import annotations`).
"""

import unittest

from neighbor_drill import (
    GENETIC_CODE_3LETTER,
    translate,
    SiftScorer,
    EditOption,
    editable_a_positions,
    enumerate_neighbor_edits,
    codon_cds_start,
    build_edit_hgvs,
    set_score,
    evaluate_options,
)


def _missense(positions, codon, aa, sift, baseline=0.1):
    """Build a scored missense EditOption (SIFT) for selection tests."""
    o = EditOption(tuple(positions), codon, aa)
    set_score(o, SiftScorer(), value=sift,
              label=("tolerated" if sift >= 0.05 else "deleterious"), baseline=baseline)
    return o


class TestEnumeration(unittest.TestCase):

    def test_editable_a_positions(self):
        self.assertEqual(editable_a_positions("CGU"), [])          # no adenosine
        self.assertEqual(editable_a_positions("CAU"), [2])
        self.assertEqual(editable_a_positions("CAA"), [2, 3])
        self.assertEqual(editable_a_positions("AAA"), [1, 2, 3])
        # accepts DNA + lowercase
        self.assertEqual(editable_a_positions("caa"), [2, 3])
        self.assertEqual(editable_a_positions("CAT"), [2])

    def test_option_counts(self):
        self.assertEqual(len(enumerate_neighbor_edits("CGA", "Arg")), 1)   # 1 A -> 1 option
        self.assertEqual(len(enumerate_neighbor_edits("CAA", "Gln")), 3)   # 2 A -> 3 options
        self.assertEqual(len(enumerate_neighbor_edits("AAA", "Lys")), 7)   # 3 A -> 7 options
        self.assertEqual(enumerate_neighbor_edits("CGU", "Arg"), [])       # 0 A -> nothing

    def test_ca_example_from_plan(self):
        # C>A: ref CAA (Gln) -> alt AAA (Lys). 3 editable A's -> 7 options, none restore Gln.
        opts = enumerate_neighbor_edits("AAA", ref_aa="Gln")
        by_combo = {o.edited_positions: o for o in opts}
        self.assertEqual((by_combo[(1,)].post_codon, by_combo[(1,)].post_aa), ("GAA", "Glu"))
        self.assertEqual((by_combo[(2,)].post_codon, by_combo[(2,)].post_aa), ("AGA", "Arg"))
        self.assertEqual((by_combo[(3,)].post_codon, by_combo[(3,)].post_aa), ("AAG", "Lys"))
        self.assertEqual((by_combo[(1, 2, 3)].post_codon, by_combo[(1, 2, 3)].post_aa), ("GGG", "Gly"))
        self.assertTrue(all(o.restores_reference is False for o in opts))

    def test_no_edit_can_create_a_stop(self):
        # An A->G edit of a sense codon must never yield a stop ('Ter').
        for codon, aa in GENETIC_CODE_3LETTER.items():
            if aa == "Ter":
                continue  # start only from sense codons
            for o in enumerate_neighbor_edits(codon, aa):
                self.assertNotEqual(o.post_aa, "Ter", f"{codon} -> {o.post_codon} became a stop")

    def test_restores_reference_flag(self):
        # Mechanism check: option {3} = AGG = Arg should be flagged when ref_aa == 'Arg'.
        opts = enumerate_neighbor_edits("AGA", ref_aa="Arg")
        by_combo = {o.edited_positions: o for o in opts}
        self.assertEqual(by_combo[(3,)].post_aa, "Arg")
        self.assertTrue(by_combo[(3,)].restores_reference)
        self.assertFalse(by_combo[(1,)].restores_reference)
        self.assertFalse(by_combo[(1, 3)].restores_reference)


class TestHgvs(unittest.TestCase):

    def test_codon_cds_start(self):
        self.assertEqual(codon_cds_start(101, 2), 100)
        self.assertEqual(codon_cds_start(100, 1), 100)
        self.assertEqual(codon_cds_start(102, 3), 100)

    def test_build_edit_hgvs_single_sub(self):
        # ref CAU, post CGU: one base differs at codon pos 2 -> c.101A>G (DNA alphabet)
        h = build_edit_hgvs("CAU", "CGU", cds_start=100, enst_ver="ENST00000123456.7")
        self.assertEqual(h, "ENST00000123456.7:c.101A>G")

    def test_build_edit_hgvs_delins(self):
        # two+ bases differ -> delins over the whole codon, emitted in DNA (T)
        h = build_edit_hgvs("CAU", "CGG", cds_start=100, enst_ver="ENST1.1")
        self.assertEqual(h, "ENST1.1:c.100_102delinsCGG")


class TestSelection(unittest.TestCase):

    def test_restore_beats_better_sift(self):
        # A restorer (no SIFT) must win is_best over a tolerated, high-SIFT missense option.
        restorer = EditOption((3,), "AAG", "Lys", restores_reference=True)
        great = _missense([2], "CGU", "Arg", sift=0.9)
        evaluate_options([great, restorer], SiftScorer(), variant_position=1)
        self.assertTrue(restorer.is_best)
        self.assertFalse(great.is_best)
        self.assertTrue(restorer.improves)        # both 'improve'; only one is best
        self.assertTrue(great.improves)

    def test_multiple_restorers(self):
        r1 = EditOption((2,), "AAG", "Lys", restores_reference=True)      # 1 edit
        r2 = EditOption((2, 3), "AGG", "Lys", restores_reference=True)    # 2 edits
        miss = _missense([1], "GAA", "Glu", sift=0.8)
        evaluate_options([r2, miss, r1], SiftScorer(), variant_position=1)
        self.assertTrue(r1.improves and r2.improves and miss.improves)   # all flagged
        self.assertTrue(r1.is_best)                                      # fewest edits among restorers
        self.assertFalse(r2.is_best)
        self.assertFalse(miss.is_best)

    def test_tiebreak_fewest_then_variant_position(self):
        # Same SIFT for all; variant_position = 1.
        var1 = _missense([1], "GAA", "Glu", sift=0.7)        # 1 edit, edits the variant position
        nbr1 = _missense([2], "AGA", "Arg", sift=0.7)        # 1 edit, neighbour
        nbr2 = _missense([2, 3], "AGG", "Arg", sift=0.7)     # 2 edits
        evaluate_options([nbr2, nbr1, var1], SiftScorer(), variant_position=1)
        self.assertTrue(var1.is_best)                        # fewest-edits tie -> variant wins
        self.assertFalse(nbr1.is_best)
        self.assertFalse(nbr2.is_best)

    def test_tiebreak_fewest_edits_overrides_variant(self):
        # A single-neighbour edit beats a (variant+neighbour) 2-edit option at equal score.
        var2 = _missense([1, 2], "GGA", "Gly", sift=0.7)     # 2 edits, includes variant
        nbr1 = _missense([3], "AAG", "Lys", sift=0.7)        # 1 edit, neighbour
        evaluate_options([var2, nbr1], SiftScorer(), variant_position=1)
        self.assertTrue(nbr1.is_best)
        self.assertFalse(var2.is_best)

    def test_improves_requires_tolerated_and_better(self):
        not_tolerated = _missense([2], "AGA", "Arg", sift=0.02, baseline=0.001)  # better but < 0.05
        worse = _missense([2], "AGA", "Arg", sift=0.5, baseline=0.9)             # tolerated but worse
        good = _missense([2], "AGA", "Arg", sift=0.5, baseline=0.1)              # tolerated + better
        for o in (not_tolerated, worse, good):
            evaluate_options([o], SiftScorer(), variant_position=1)
        self.assertFalse(not_tolerated.improves)
        self.assertFalse(worse.improves)
        self.assertTrue(good.improves)

    def test_no_scores_no_best(self):
        # No restorers and no scores -> nothing selected, no crash.
        o = EditOption((2,), "AGA", "Arg")
        evaluate_options([o], SiftScorer(), variant_position=1)
        self.assertFalse(o.is_best)
        self.assertFalse(o.improves)

    def test_set_score_computes_improvement(self):
        o = EditOption((2,), "AGA", "Arg")
        set_score(o, SiftScorer(), value=0.8, label="tolerated", baseline=0.1)
        self.assertEqual(o.scores["SIFT"].improvement, 0.7)
        # missing baseline or value -> no improvement
        o2 = EditOption((2,), "AGA", "Arg")
        set_score(o2, SiftScorer(), value=None, label=None, baseline=0.1)
        self.assertIsNone(o2.scores["SIFT"].value)
        self.assertIsNone(o2.scores["SIFT"].improvement)


class TestTranslate(unittest.TestCase):

    def test_translate_accepts_dna_and_rna(self):
        self.assertEqual(translate("AUG"), "Met")
        self.assertEqual(translate("ATG"), "Met")
        self.assertEqual(translate("atg"), "Met")
        self.assertEqual(translate("UAA"), "Ter")
        self.assertEqual(translate("XYZ"), "Unknown")


if __name__ == "__main__":
    unittest.main()
