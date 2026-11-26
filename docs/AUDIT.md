# What the two audits do

Both scripts are adversarial: they were written to catch the author, not to confirm him. They
read the deposited manuscript and the deposited tables, and they fail loudly when a number in
the text cannot be found in a table or disagrees after rounding.

## 345 — text and tables versus data

126 checks. Each is a triple (location in the manuscript, expected value, source table and
column). The script extracts the value from the DOCX, recomputes it from the deposited table, and
compares within the rounding used in the text.

Coverage: Table 1 (rho_obs, rho_corr, intervals), Table 2 ranges, the readability counts at every
threshold in Table S3, the group contrast and its difference-in-differences under both denominators
(12.266 fixed, 12.265 covered), the sample-level bootstrap interval (9.4-14.6), the 39 assessable
samples and the 10 readable ones, the 29 with r_full,cov 0.003-0.51, the coverage-aware reliability
over the 41 samples, the double-zero decomposition (5,257 copies, 27.4% covered on both haplotypes,
2.0-92% per sample), the phasing diagnostics (Hap1 at 88.3% with neither haplotype covered,
GM12878 hap2 at 28.0%), the six pluripotent pairs (0.92-0.98 against H1's 0.93-0.98), the
allele-level correlations, and the group sizes the text quotes (26 lymphoblastoid lines, 25 of them
assessable).

A few checks exist only to stop a phrase from coming back: that no sentence attributes the c5
difference to sequencing depth (no depth analysis was run), that the double-zero sentence keeps its
direction, and that the 1,032-copy panel is not given a copy-level interval.

The last revision added format guards rather than number guards: every supplementary citation must
carry its Online Resource number, no figure caption may end with a punctuation mark, the compendium
identifiers (Hek293T, Panc1, Hap1) must not reappear, the seven abbreviations must be written out at
first use, the drawing software and its version must be stated in Methods, and the eight declaration
headings must all be present.

Last run against the deposited manuscript: **126 checks, 0 mismatch**.

## 346 — figures versus text

97 checks. The values a figure draws are written by the plotting script to
`figures/source_data/fig*_numbers.json`; this script matches every one of them against the text.
It is the check that catches a panel and a sentence drifting apart after an edit.

Last run against the deposited manuscript: **97 checks, 0 mismatch**.

## Known limitation

Neither script can verify the upstream extraction step (bigWig to per-copy pa5), because that
needs the raw tracks. The independent guard for that step is inside `340_fig_common.py`:
`spot_check()` recomputes pa5 for three random copies from the public bigWig and compares it with
the deposited table at every figure build.

pa5 is a mean over a fixed 1,000-bp window, so a window with no fibre coverage is scored as zero
in the same way as a window measured as closed. The coverage tracks separate the two cases, and
`t15_table_s4_cov.tsv` reports the split-half correlation obtained when a copy is kept on coverage
(both haplotypes >= 800 bp) rather than on a fixed denominator. Window-level coverage in the pooled
tracks is close to complete (0.02% of 42,312 window-by-sample pairs leave a base uncovered), so the
two denominators give the same difference-in-differences to within 0.01 points; the difference
between them is in the haplotype tracks, where 9.1% of copy-haplotypes carry no coverage at all.
