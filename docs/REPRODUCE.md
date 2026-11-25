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
python 333_t15_tables.py      # Tables S3, S5, S7 (+ the atlas behind S4)
python 337_t15_s8.py          # Table S8
python 356_online_resource.py # assembles Tables S1-S8 into one xlsx
```

`356` rebuilds `manuscript/Online_Resource_1_Tables_S1-S8.xlsx` cell-for-cell.

## 3. Audits

```
python 345_text_audit.py      # 70 checks: text and tables vs data/derived
python 346_figure_audit.py    # 87 checks: *_numbers.json vs text
```

Both write TSVs to `docs/audit/` and print a mismatch count. Expected: 0 and 0.

## 4. Reference metadata

```
python 359_doi_check.py
```

Queries Crossref for every reference and prints title/volume/page agreement. Four references
report WARN for metadata style rather than content (journal abbreviation, first page only, a
preprint with no journal, and one DOI whose Crossref record carries no page range but whose
Europe PMC record does). See `docs/AUDIT.md`.

## 5. What was checked, and what could not be

* Verified: figure values, table values, the numbers quoted in the text, the OR1 workbook,
  reference metadata.
* Not verified here: the upstream extraction from bigWig (needs the ~hundreds of GB of public
  tracks), and the wet-lab provenance of the compendium itself, which belongs to the original
  publications.
