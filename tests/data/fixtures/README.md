# Test fixtures

Real variants and their **real** VEP output, sliced from the paper run.
Regenerate with `python scripts/build_fixtures.py`; verify with `--check`.

Selection is by criterion (lowest matching variant id), never hand-picked, so this set
is reproducible and each entry documents why it is here.

`vep_slice.txt` keeps the original `##` header, so the VEP version and the bundled
ClinVar/SIFT/dbSNP releases travel with the fixture (see `docs/provenance.md`).

| variant (chr:pos:ref:alt) | criterion | why it is in the fixture |
|---|---|---|
| `11:4566422:T:C` | start_lost | rare consequence class (n=10) |
| `1:103616:CTTTTTTTT:CTTTT` | deletion | class != SNV; exercises the CADD SNV-only gate |
| `1:10725142:G:C` | multi_consequence | several VEP terms on one variant; exercises consequence collapsing and most-severe ranking |
| `1:2234752:G:A` | g2a_plus_strand + missense_alt_codon_has_A | directly correctable, plus strand: genomic G>A == coding G>A; also: the Improve path: an editable adenosine inside the mutant codon |
| `1:2234766:C:T` | stopgained + codon_position_1 | the Rescue path: stop codon recoded to TGG; also: codon offset boundary |
| `1:2237536:G:T` | improve_not_also_fix | residue optimisation with no simpler route: missense, editable adenosine in the mutant codon, and NOT directly revertible |
| `1:46658987:C:T` | missing_sift | missense with no SIFT score; must degrade, not crash |
| `1:840685:G:GT` | insertion | class != SNV; ref/alt lengths differ, end != start |
| `1:8418724:C:T` | c2t_minus_strand | directly correctable, minus strand: genomic C>T == coding G>A. The case the 218-vs-138 error got wrong |
| `1:8418972:G:A` | g2a_on_minus_not_correctable | NOT correctable despite being genomic G>A: on a minus-strand gene it reads as C>T. The complement of the case above; both are needed to pin the strand rule |
| `1:8421289:G:A` | clinvar_pathogenic | ClinVar pathogenic under the Figure 1C rule |
| `1:865654:C:T` | synonymous + codon_position_3 | in the CDS set but not treatable; also: codon offset boundary |
| `1:878673:A:T` | no_rs_id | no Existing_variation; rs_id must be NULL not empty |
| `1:878759:T:G` | splice | splice-site consequence bucket |
| `X:601587:G:A` | chrX | non-autosome; guards chromosome handling |
