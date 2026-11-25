# What the two audits do

Both scripts are adversarial: they were written to catch the author, not to confirm him. They
read the deposited manuscript and the deposited tables, and they fail loudly when a number in
the text cannot be found in a table or disagrees after rounding.

## 345 — text and tables versus data

Seventy-one checks. Each is a triple (location in the manuscript, expected value, source table
and column). The script extracts the value from the DOCX, recomputes it from the deposited table,
and compares within the rounding used in the text.

Coverage: Table 1 (rho_obs, rho_corr and intervals), Table 2 ranges, the stability of the readable
count between 0.50 and 0.667 (Table S3 holds all seven thresholds), the group contrasts
(+0.71 / +0.59 / +0.90), the 14 cells with r_s_nz <= 0, the phasing diagnostics (including
GM12878 hap2 at 28.0%), the allele-level correlations, and the per-sample values quoted in the
Results.

Last run against the deposited manuscript: **71 checks, 0 mismatch**.

## 346 — figures versus text

Eighty-seven checks. The values a figure draws are written by the plotting script to
`figures/source_data/fig*_numbers.json`; this script matches every one of them against the text.
It is the check that catches a panel and a sentence drifting apart after an edit.

Last run against the deposited manuscript: **87 checks, 0 mismatch**.

## Known limitation

Neither script can verify the upstream extraction step (bigWig to per-copy pa5), because that
needs the raw tracks. The independent guard for that step is inside `340_fig_common.py`:
`spot_check()` recomputes pa5 for three random copies from the public bigWig and compares it with
the deposited table at every figure build.

pa5 is a mean over a fixed 1,000-bp window, so a window with no fibre coverage is scored as zero
in the same way as a window measured as closed. Nothing in the deposited tables separates those
two cases; the coverage tracks do, and that recomputation is the one piece of work listed under
"Open items" below.

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

## Open items

The coverage-aware recomputation of split-half reliability (`r_s_cov`, copies retained when both
haplotype tracks have covered bases) is not yet deposited. It can be computed for the 25 samples
that carry a coverage track; for the four samples that moved from readable to not readable
(Panc1, A549, THP-1, HG002) no coverage track exists, so the manuscript reports only the bound
that up to three of them would return to readable if their double-zero copies are closed rather
than unphased.
