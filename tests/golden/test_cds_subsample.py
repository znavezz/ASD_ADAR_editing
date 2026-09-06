"""Golden tests for the committed CDS subsample.

`tests/data/golden/cds_subsample.csv` is a 200-row slice of the 9,962-row snapshot that
every figure after Figure 1 is built from. The full snapshot lives under the gitignored
`Output/` tree, so this subsample is the committed, value-level regression fixture: a
refactor that changes individual annotations -- not just row counts -- fails here with a
readable diff.

These tests assert intrinsic properties of the committed file, so they need no database,
no network, and no access to the full snapshot.
"""
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SUBSAMPLE = REPO_ROOT / "tests" / "data" / "golden" / "cds_subsample.csv"

#: Column set and order of the snapshot, locked. A change here means the query behind
#: `Output/Intermediate_tables/00_db_snapshot/` changed shape and the figures may move.
EXPECTED_COLUMNS = [
    "id", "REF.VEP", "ALT.VEP", "class", "cds_position", "protein_position", "exon",
    "intron", "inheritance", "rs_id", "nmd_escaping_variant", "position_in_codon", "chr",
    "start", "end", "STRAND", "gene_ensg", "gene_symbol", "gene_strand", "ref_codon_nt1",
    "ref_codon_nt2", "ref_codon_nt3", "ref_codon_stop", "ref_aa", "ref_aa_name",
    "alt_codon_nt1", "alt_codon_nt2", "alt_codon_nt3", "alt_codon_stop", "alt_aa",
    "alt_aa_name", "Consequence", "New_Consequence", "patient_count", "paper_keys",
    "CLIN_SIG", "SIFT_pred", "am_class", "condel_pred", "polyphen_pred", "SIFT_score",
    "am_pathogenicity", "polyphen_score", "REVEL", "BLOSUM62", "condel_score", "pLI",
    "PhastCons100", "LoFtool", "GENE_PHENO", "CADD_phred", "sfari_syndromic",
    "sfari_gene_score", "sfari_reports", "REF", "ALT",
]

#: The five consequence classes in the snapshot manifest.
EXPECTED_CONSEQUENCES = {"Missense", "Nonsense", "Splice", "Start_lost", "Synonymous"}


@pytest.fixture(scope="module")
def sample() -> pd.DataFrame:
    if not SUBSAMPLE.exists():
        pytest.fail(f"{SUBSAMPLE} is missing; regenerate with `python scripts/cut_subsample.py`.")
    return pd.read_csv(SUBSAMPLE, low_memory=False)


class TestShape:
    def test_row_count(self, sample):
        assert len(sample) == 200

    def test_columns_are_locked(self, sample):
        assert list(sample.columns) == EXPECTED_COLUMNS

    def test_ids_are_unique_and_sorted(self, sample):
        ids = sample["id"].tolist()
        assert len(set(ids)) == len(ids)
        assert ids == sorted(ids)

    def test_stays_small_enough_to_commit(self):
        assert SUBSAMPLE.stat().st_size < 512 * 1024


class TestCoverage:
    """The sample must exercise every shape a refactor could plausibly break."""

    def test_every_consequence_class_present(self, sample):
        assert EXPECTED_CONSEQUENCES <= set(sample["New_Consequence"].dropna())

    def test_both_gene_strands_present(self, sample):
        assert set(sample["gene_strand"].dropna()) == {"+", "-"}

    def test_all_three_codon_positions_present(self, sample):
        assert set(sample["position_in_codon"].dropna().astype(int)) == {1, 2, 3}

    def test_adar_fix_patterns_present(self, sample):
        """G>A on the plus strand and C>T on the minus strand -- the two variant shapes
        that A-to-I editing reverts. Every editability classification keys off these."""
        plus = (sample["REF"] == "G") & (sample["ALT"] == "A") & (sample["gene_strand"] == "+")
        minus = (sample["REF"] == "C") & (sample["ALT"] == "T") & (sample["gene_strand"] == "-")
        assert plus.any(), "no G>A-on-plus-strand row"
        assert minus.any(), "no C>T-on-minus-strand row"
        assert (~(plus | minus)).any(), "no non-Fix row"

    def test_missing_and_present_sift_both_covered(self, sample):
        assert sample["SIFT_score"].isna().any()
        assert sample["SIFT_score"].notna().any()

    def test_multi_consequence_row_present(self, sample):
        assert sample["Consequence"].str.contains(",", na=False).any()

    def test_nmd_escaping_both_values_present(self, sample):
        assert set(sample["nmd_escaping_variant"].astype(bool)) == {True, False}

    def test_intron_and_missing_cds_position_covered(self, sample):
        assert sample["intron"].notna().any()
        assert sample["cds_position"].isna().any()

    def test_every_sfari_gene_score_present(self, sample):
        assert set(sample["sfari_gene_score"].dropna().astype(int)) == {1, 2, 3}


class TestInvariants:
    """Properties that must hold for every row, in the sample and in the full table."""

    def test_coordinates_are_one_based_and_ordered(self, sample):
        """Positions are 1-based inclusive, so start >= 1 and end >= start."""
        assert (sample["start"] >= 1).all()
        assert (sample["end"] >= sample["start"]).all()

    def test_snv_end_equals_start(self, sample):
        snv = sample[sample["class"] == "SNV"]
        assert not snv.empty
        assert (snv["end"] == snv["start"]).all()

    def test_ref_and_alt_differ(self, sample):
        assert (sample["REF"] != sample["ALT"]).all()

    def test_alleles_are_unambiguous(self, sample):
        assert sample["REF"].str.fullmatch(r"[ACGT]+").all()
        assert sample["ALT"].str.fullmatch(r"[ACGT]+").all()

    def test_position_in_codon_is_valid(self, sample):
        present = sample["position_in_codon"].dropna()
        assert present.isin([1, 2, 3]).all()

    def test_sift_score_is_a_probability(self, sample):
        present = sample["SIFT_score"].dropna()
        assert ((present >= 0) & (present <= 1)).all()

    def test_codon_columns_are_single_nucleotides(self, sample):
        for column in ["ref_codon_nt1", "ref_codon_nt2", "ref_codon_nt3",
                       "alt_codon_nt1", "alt_codon_nt2", "alt_codon_nt3"]:
            present = sample[column].dropna()
            assert present.isin(list("ACGT")).all(), f"{column} has non-ACGT values"

    def test_cds_position_agrees_with_codon_position(self, sample):
        """VEP's CDS coordinate and position-in-codon must be consistent: for a 1-based
        CDS position p, the offset within its codon is ((p - 1) mod 3) + 1."""
        rows = sample.dropna(subset=["cds_position", "position_in_codon"])
        assert not rows.empty
        derived = ((rows["cds_position"].astype(int) - 1) % 3) + 1
        assert (derived == rows["position_in_codon"].astype(int)).all()
