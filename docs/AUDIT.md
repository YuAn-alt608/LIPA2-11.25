# What the two audits do

Both scripts are adversarial: they were written to catch the author, not to confirm him. They
read the deposited manuscript and the deposited tables, and they fail loudly when a number in
the text cannot be found in a table or disagrees after rounding.

## 345 — text and tables versus data

Seventy checks. Each is a triple (location in the manuscript, expected value, source table and
column). The script extracts the value from the DOCX, recomputes it from the deposited table,
and compares within the rounding used in the text.

Coverage: Table 1 (rho_obs, rho_corr and intervals), Table 2 ranges, the threshold counts
(10/10/10/8/5/4/2), the group contrasts (+0.71 / +0.59 / +0.90), the 14 cells with r_s_nz <= 0,
the phasing diagnostics (including GM12878 hap2 at 28.0%), the allele-level correlations, and the
per-sample values quoted in the Results.

Last run against the deposited manuscript: **70 checks, 0 mismatch**.

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
