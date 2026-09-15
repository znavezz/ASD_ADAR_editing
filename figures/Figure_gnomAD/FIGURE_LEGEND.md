# Figure legend

Text for the manuscript. Everything here was deliberately kept off the plot itself, so the
panel carries only the data and the legend carries the qualifications.

---

**Figure X. Population allele frequency of variants amenable to ADAR-mediated RNA editing.**
Distribution of gnomAD exome allele frequency across the three editing classes: directly
correctable G>A variants (Direct Repair, n = 2,421), non-G>A nonsense variants amenable to
stop-codon recoding (Nonsense Rescue, n = 894), and non-G>A missense variants amenable to
alternative codon editing (Missense Optimization, n = 158). Bars show the proportion of each
class falling in each frequency band; bands are ordered from rarest (left) to most common
(right). Variants not observed in gnomAD are shown separately in grey rather than as a
frequency, since absence from a population reference is evidence of rarity rather than a
missing measurement. Allele frequencies are gnomAD exomes r2.1, obtained from the Ensembl VEP
GRCh37 cache, which does not provide gnomAD genome frequencies. Across all three classes the
large majority of variants are absent from gnomAD or present below 0.01%: 91.6% for Direct
Repair, 97.2% for Nonsense Rescue and 92.4% for Missense Optimization. Seventeen Direct Repair
variants (0.7%) exceed a frequency of 1%, of which two exceed 5% - the threshold at which
population frequency alone is considered evidence against pathogenicity (ACMG BA1). Neither
amino acid-level class contains a variant above 1%.

---

## Numbers quoted above, for checking

| Class | n | Absent | < 0.0001 | 0.0001-0.001 | 0.001-0.01 | >= 0.01 | absent or <1e-4 |
|---|---|---|---|---|---|---|---|
| Direct Repair | 2,421 | 1,058 | 1,160 | 144 | 42 | 17 | **91.6%** |
| Nonsense Rescue | 894 | 669 | 200 | 20 | 5 | 0 | **97.2%** |
| Missense Optimization | 158 | 93 | 53 | 8 | 4 | 0 | **92.4%** |

Of the 17 Direct Repair variants at >= 0.01, two are at >= 0.05.

Regenerate with `Rscript figures/Figure_gnomAD/00_af_distribution.R`; the per-band counts are
written alongside the figure as `Figure_gnomAD_AF_distribution_data.csv`.

## If this becomes a manuscript figure

Two sentences would belong in the Methods, since neither is visible from the panel:

  * Allele frequencies were taken from the gnomAD exome fields of the Ensembl VEP annotation
    (r2.1, GRCh37 cache). Genome frequencies are not available in that cache.
  * Variants absent from gnomAD were retained and analysed as a distinct category; allele
    frequency was not used to filter any variant set.
