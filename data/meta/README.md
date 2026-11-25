# Sample manifest

`fiberseq_samples.tsv` lists the 41 samples analysed in the paper:

| column | meaning |
| --- | --- |
| `sample` | sample name used throughout the manuscript |
| `ps_id` | compendium identifier |
| `tissue` | cell type or tissue as annotated by the compendium |
| `default` | whether the compendium marks the sample as a default track |
| `percent_accessible_all` | hg38 bigWig, all fibres |
| `percent_accessible_hap1` | hg38 bigWig, haplotype 1 |
| `percent_accessible_hap2` | hg38 bigWig, haplotype 2 |
| `peaks` | FIRE peak bigBed |

The URLs point at the public UCSC/PacBio Fiber-seq Compendium. They are the exact files the
analysis read; `340_fig_common.py` uses them for the per-build spot check.
