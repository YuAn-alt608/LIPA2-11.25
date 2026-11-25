# What the two audits do

Both scripts are adversarial: they were written to catch the author, not to confirm him. They
read the deposited manuscript and the deposited tables, and they fail loudly when a number in
the text cannot be found in a table or disagrees after rounding.

## 345 — text and tables versus data

Eighty-five checks. Each is a triple (location in the manuscript, expected value, source table
and column). The script extracts the value from the DOCX, recomputes it from the deposited table,
and compares within the rounding used in the text.

Coverage: Table 1 (rho_obs, rho_corr and intervals), Table 2 ranges, the stability of the readable
count between 0.50 and 0.667 (Table S3 holds all seven thresholds), the group contrasts
(+0.71 / +0.59 / +0.90), the 14 cells with r_s_nz <= 0, the phasing diagnostics (including
GM12878 hap2 at 28.0%), the allele-level correlations, and the per-sample values quoted in the
Results, and the coverage-aware reliability figures added in r27 (41 samples with a coverage
track, 39 assessable, 27.8% of removed copies covered on both haplotypes, r_s_cov 0.001-0.34 over
the 29 not readable samples).

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
rather than removed for being double-zero. All 41 samples carry such a track; for Hap1 the split
half is not defined and GM12878 lacks a usable second haplotype, so the 39 assessable samples are
the ones on which the split is interpreted.

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

## Coverage-aware reliability (deposited in r26, closed in r27)

`t15_table_s4_cov.tsv` adds four columns to the atlas: the number of copies covered on both
haplotypes, the fraction of double-zero copies that carry such coverage, and the split-half
correlation `r_s_cov` with its Spearman-Brown value.

r26 could compute these columns for 26 samples, because the 15 samples fetched earliest had been
entered through a path that did not record `covered_bp`. r27 fills them in: `381_fetch_cov15.py`
re-reads both haplotype tracks of those 15 samples and keeps only the coverage column, the pa5
column staying as deposited. Each fetched track is checked against the deposited pa5 before it is
accepted, and all 30 tracks agree to 1e-15, so the coverage values are aligned window-for-window
with the numbers already in the atlas.

With the full set, 5,257 double-zero copies carry coverage, and 1,463 of them (27.8%) are covered
on both haplotypes. `r_s_cov` runs from 0.001 to 0.34 over the 29 not readable samples and from
0.60 to 0.89 over the ten readable ones, so no sample crosses the threshold and the readable
classification is unchanged. Within the 25 lymphoblastoid lines it stays at 0.001-0.14. None of
the 39 assessable samples has `r_s_cov` at or below zero, which is the evidence that the negative
values in the double-zero column belong to the estimator.

The four samples that moved from readable to not readable are the case the r26 text could only
bound. Only 2-9% of their double-zero copies are covered on both haplotypes, so those copies lack
reads rather than being measured closed, and all four stay below the threshold with
`r_full_cov` 0.22-0.51. The bound the r26 text gave ("up to three would return") is therefore
resolved downward: none returns.
