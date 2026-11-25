# What the two audits do

Both scripts are adversarial: they were written to catch the author, not to confirm him. They
read the deposited manuscript and the deposited tables, and they fail loudly when a number in
the text cannot be found in a table or disagrees after rounding.

## 345 — text and tables versus data

Eighty checks. Each is a triple (location in the manuscript, expected value, source table
and column). The script extracts the value from the DOCX, recomputes it from the deposited table,
and compares within the rounding used in the text.

Coverage: Table 1 (rho_obs, rho_corr and intervals), Table 2 ranges, the stability of the readable
count between 0.50 and 0.667 (Table S3 holds all seven thresholds), the group contrasts
(+0.71 / +0.59 / +0.90), the 14 cells with r_s_nz <= 0, the phasing diagnostics (including
GM12878 hap2 at 28.0%), the allele-level correlations, and the per-sample values quoted in the
Results, and the coverage-aware reliability figures added in r26 (26 samples with a coverage
track, 25 assessable, 56% of removed copies covered on both haplotypes).

Last run against the deposited manuscript: **80 checks, 0 mismatch**.

## 346 — figures versus text

Ninety-one checks. The values a figure draws are written by the plotting script to
`figures/source_data/fig*_numbers.json`; this script matches every one of them against the text.
It is the check that catches a panel and a sentence drifting apart after an edit.

Last run against the deposited manuscript: **91 checks, 0 mismatch**.

## Known limitation

Neither script can verify the upstream extraction step (bigWig to per-copy pa5), because that
needs the raw tracks. The independent guard for that step is inside `340_fig_common.py`:
`spot_check()` recomputes pa5 for three random copies from the public bigWig and compares it with
the deposited table at every figure build.

pa5 is a mean over a fixed 1,000-bp window, so a window with no fibre coverage is scored as zero
in the same way as a window measured as closed. The coverage tracks separate the two cases, and
`t15_table_s4_cov.tsv` reports the split-half correlation obtained when a copy is kept on coverage
rather than removed for being double-zero; the limitation applies to the 15 samples that carry no
coverage track, where the two cases remain indistinguishable.

## Earlier row-order fault, kept for the record

An early version of the pipeline split rows into 24 chunks (one per thread) and concatenated them
in chunk order, which permuted the row order of the per-copy tables by a fixed stride of 24.
Within-table statistics were unaffected; joins by copy identifier across files were wrong. The
permutation was recovered by matching H1 ranks, every copy was recomputed directly from the
source bigWig (1,032 L1PA2 copies and 7,222 full-length L1PA copies; maximum difference
7.1e-15 against the corrected table, 28.6 percentage points against the faulty one), and the
affected analyses were re-run. The fault recurred in a later fetch of haplotype tracks for 26
samples, where the identifier column stayed in the original order so the identifier assertions
passed; it was caught because each sample's pooled track and the mean of its two haplotype tracks
became uncorrelated, and the 26 samples were re-fetched value-by-value. The manuscript reports
only the resulting validation, because this is an engineering history rather than a property of
the measurement.

## Coverage-aware reliability (deposited in r26)

`t15_table_s4_cov.tsv` adds four columns to the atlas: the number of copies covered on both
haplotypes, the fraction of double-zero copies that carry such coverage, and the split-half
correlation `r_s_cov` with its Spearman-Brown value. For the 25 assessable samples that carry a
coverage track, `r_s_cov` runs from 0.001 to 0.209 with none at or below zero and none above 0.5,
so the readable classification does not change and the negative values in the double-zero column
are a property of the estimator rather than of the tracks. The four samples that moved from
readable to not readable (Panc1, A549, THP-1, HG002) carry no coverage track, so for them the
manuscript reports only the bound: up to three would return to readable if their double-zero
copies are closed rather than unphased.
