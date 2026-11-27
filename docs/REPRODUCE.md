# Reproducing the paper

Everything below was executed on a fresh copy of this repository before deposit. Python 3.11.2,
the versions in `requirements-lock.txt`, Windows 11. Commands are run from `code/`.

## 1. Panels

```
python 341_fig1.py     # Fig. 1, 174 mm wide
python 342_fig2.py     # Fig. 2
python 343_fig3.py     # Fig. 3
python 344_figS.py     # Fig. S1-S5 in one pass
```

Output goes to `figures/generated/` (PDF, EPS, PNG at ~600 dpi for 174 mm, plus
`*_numbers.json`). The deposited panels in `figures/` were produced by exactly these scripts;
`figures/generated/` is empty in the clone so that a rebuild cannot be confused with the deposit.

The eight `*_numbers.json` files are the bridge between code and text: every value the figures
draw is written there, and `346_figure_audit.py` matches those values against the manuscript.

Offline: `LIPA2_SKIP_SPOTCHECK=1` skips the remote bigWig spot check.

## 2. Tables

```
python 387_fetch_pool_cov.py  # pooled coverage per window, 41 samples
python 389_reliability_cov.py # the coverage-aware reliability of all 41 samples
python 390_did_cov.py         # the two denominators of the difference-in-differences
python 391_table_s4_cov.py    # Table S4, coverage-aware columns (+ 800 bp floor)
python 399_t15_s5_cov.py      # Table S5, pairwise agreement corrected with r_full,cov
python 397_sim_band_cov.py    # the noise-only band on r_full,cov
python 398_threshold_cov.py   # readability under each r_full,cov threshold (Table S3)
python 333_t15_tables.py      # Tables S3, S5, S7 (+ the atlas behind S4)
python 381_fetch_cov15.py     # haplotype coverage tracks for the 15 early samples
python 375_table_s4_cov.py    # the earlier Table S4 columns (kept as *_dz)
python 337_t15_s8.py          # Table S8
python 356_online_resource.py # assembles Tables S1-S8 into one xlsx
```

`356` rebuilds `manuscript/Online_Resource_1_Tables_S1-S8.xlsx` cell-for-cell. `387` and `381`
read the public bigWig over HTTP and need network access; their output is deposited, so
everything after them runs offline.

## 3. Audits

```
python 345_text_audit.py      # 130 checks: text and tables vs data/derived
python 346_figure_audit.py    # 97 checks: *_numbers.json vs text
```

The first takes the manuscript path as an argument (default `manuscript/`). Both write TSVs to
`docs/audit/` and print a mismatch count. Expected: 0 and 0.

## 4. Reference metadata

```
python 359_doi_check.py       # every reference resolved through Crossref
```

## 5. What this cannot reproduce

The manuscript text itself. The DOCX was typeset from the deposited tables and audited against
them, but the prose was written and edited by hand; no script regenerates it. The typesetting
scripts perform no analysis and are not in this repository.
