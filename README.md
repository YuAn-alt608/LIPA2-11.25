# Reliability-aware Fiber-seq readout of young LINE-1 (L1PA2) promoter accessibility in human pluripotent cells

Frozen analysis snapshot for a manuscript submitted to *Functional & Integrative Genomics*
(article type: Methodology), by Zizheng Yuan (University of Lisbon).

All input data are public and no new sequencing was performed. This repository contains the
manuscript, the eight figure panels at print resolution, every derived table behind the reported
numbers, and the code that regenerates them.

## Reproducibility, as verified

Run from a fresh clone, using only the deposited tables:

| Check | Command | Result |
| --- | --- | --- |
| Eight figure panels | `python code/341_fig1.py` … `344_figS.py` | all eight `*_numbers.json` rewritten **byte-identically** |
| Online Resource 1 (Tables S1–S8) | `python code/356_online_resource.py` | rebuilt **cell-for-cell identical** (9 sheets) |
| Every number in the text | `python code/345_text_audit.py` | 112 checks, **0 mismatch** |
| Every number on the figures | `python code/346_figure_audit.py` | 97 checks, **0 mismatch** |

The audits read the deposited manuscript and the deposited tables; nothing is hand-entered.
See `docs/AUDIT.md` for what each check compares.

## Contents

```
manuscript/   main text (docx + pdf), cover letter, Online Resource 1 (Tables S1-S8)
figures/      Fig1-Fig3 and FigS1-FigS5 as PDF (vector), EPS and 600 dpi PNG
              source_data/  per-panel values and the *_numbers.json read by the figure code
data/         derived/  the frozen analysis tables (see docs/DATA_DICTIONARY.md)
              meta/     the 41-sample manifest with the public track URLs
code/         analysis, figure, table and audit scripts (see docs/SCRIPT_MAP.md)
docs/         data dictionary, script map, reproduction steps, audit description
SHA256SUMS.txt  checksum of every file in this repository
```

## Data sources

| Input | Source | Accessed through |
| --- | --- | --- |
| Fiber-seq `percent.accessible` (all, hap1, hap2) and FIRE peaks, 41 samples, hg38 | UCSC/PacBio Fiber-seq Compendium, the resource described by Vollger et al. (2026) *Cell* 189, 6613-6629 | per-sample URLs in `data/meta/fiberseq_samples.tsv` |
| hiPS6 RNA-seq and CUT&RUN | Adami et al. (2025) *Cell Genomics* 5, 100979 | GEO/ENA deposit of that study |
| Bisulfite data for the independent methylation comparison | ENCODE | ENCFF423PGW, ENCFF716NKX |
| Repeat annotation | RepeatMasker (Smit et al. 2015) | hg38 rmsk |

`data/meta/fiberseq_samples.tsv` maps each of the 41 samples to its four public track URLs, so
every per-copy value in `data/derived/` can be traced back to a byte range of a public file.

## Quick start

```
python 3.11
pip install -r requirements-lock.txt

cd code
python 341_fig1.py            # writes figures/generated/Fig1.pdf/.png and fig1_numbers.json
python 342_fig2.py
python 343_fig3.py
python 344_figS.py            # all five supplementary panels
python 356_online_resource.py # rebuilds manuscript/Online_Resource_1_Tables_S1-S8.xlsx
python 345_text_audit.py      # text vs tables -> docs/audit/
python 346_figure_audit.py    # figure values vs text -> docs/audit/
```

Each figure script starts with a spot check that recomputes pa5 for three random L1PA2 copies
from the public bigWig over HTTP and compares them with the deposited table. That check needs
network access; set `LIPA2_SKIP_SPOTCHECK=1` to skip it when offline (the figures themselves are
drawn from the deposited tables and do not need the network).

## What is not in this repository

* **Raw tracks.** The 41 samples amount to hundreds of GB of bigWig. They are public and
  addressable by the URLs in `data/meta/fiberseq_samples.tsv`, so they are not mirrored here.
* **Upstream per-copy extraction.** The step that turns bigWig into per-copy pa5 (`236_pa5_expand.py`)
  is included for provenance, but it needs the tracks. Its output, `data/derived/t9b/` and
  `data/derived/t15/t15_pa5_hap_L1PA2_41.tsv.gz`, is deposited, and everything downstream of it runs
  from the deposit alone.
* **Manuscript typesetting.** The python-docx scripts that applied journal formatting and tracked
  changes perform no analysis and are not included; the typeset manuscript is deposited as the
  output of that step.

## Licence

Code: MIT (`LICENSE`). Derived tables, figures and manuscript: CC BY 4.0 (`LICENSE-DATA`).

## Contact

Zizheng Yuan — yuan-zizheng@edu.ulisboa.pt — https://github.com/YuAn-alt608/LIPA2-11.25
