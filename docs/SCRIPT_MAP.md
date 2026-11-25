# Script map

Scripts are numbered as they were developed; gaps are steps that were superseded and removed.
All paths inside the copies are relative to this repository, so the scripts run from `code/`
without editing.

## Figures

| Script | Produces | Reads |
| --- | --- | --- |
| `340_fig_common.py` | shared style, panel layout, colour code, `spot_check()` | `data/derived/l1pa5_copy_tracks_v2.tsv.gz`, `data/derived/t9b/`, `data/meta/fiberseq_samples.tsv` |
| `341_fig1.py` | Fig. 1 (a-f) | `t1_did_perm4_c5.tsv`, `t1_did_perm4_pa5.tsv`, `t13/t13_group_contrast_sample_first.tsv`, `t13/t13_per_sample_did.tsv`, `t1d_pa5_by_copy_9.tsv.gz` |
| `342_fig2.py` | Fig. 2 (a-b) | `t15/t15_pa5_vs_rfull.tsv`, `t15/t15_rho_vs_h1_41.tsv`, `t15/t15_table_s4_cov.tsv`, `t10_sim_band.tsv` |
| `343_fig3.py` | Fig. 3 (a-e) | `t9m/t9m_partial_rho.tsv`, `t13/t13_partial_rho_allcopies.tsv`, `t9d/t9d_coupling.tsv`, `t15/`, `t10_sim_band.tsv`, `t9m/t9m_rna_by_copy.tsv.gz`, `t9b/` |
| `344_figS.py` | Fig. S1-S5 | `t1_did_perm4_c5.tsv`, `t10_sim_band.tsv`, `t15/t15_selfconsistency.tsv`, `t9i/t9i_dA_by_copy.tsv.gz`, `t10/t10_readability_map.tsv` |

## Tables

| Script | Produces | Reads |
| --- | --- | --- |
| `333_t15_tables.py` | Tables S3, S5, S7 and the reliability atlas | `t15/t15_reliability_41.tsv`, `t15/t15_pa5_hap_L1PA2_41.tsv.gz`, `t9b/` |
| `381_fetch_cov15.py` | haplotype `covered_bp` for the 15 early samples (aligned to the deposited pa5) | `data/meta/fiberseq_samples.tsv`, `t9b/t9b_pa5_hap_L1PA2_fixed.tsv.gz` |
| `375_table_s4_cov.py` | Table S4, coverage-aware columns | `t15/t15_table_s3.tsv`, `t15/t15_pa5_hap_L1PA2_41.tsv.gz`, `t15/t15_hap_cov_L1PA2_41.tsv.gz` |
| `337_t15_s8.py` | Table S8 | same, plus `l1pa5_copy_tracks_v2.tsv.gz` |
| `356_online_resource.py` | Online Resource 1 (S1-S8) | all of the above plus `sample_labels.tsv`, `t9d/`, `t9i/` |

## Analysis

| Script | What it computes |
| --- | --- |
| `236_pa5_expand.py` | per-copy pa5 from the public bigWig (upstream; needs the tracks) |
| `332_reliability_41_v2.py` | split-half and haplotype reliability for 41 samples |
| `316_reliability_41.py` | first-pass reliability with the pooled-versus-haplotype check |
| `295_dA_dM_coupling.py` | allele-level dA-dM coupling |
| `297_readability_extend.py`, `303_reliability_doublezero.py`, `304_merge_readability_nz.py` | readability rule with and without double-zero copies |
| `306_anchor_partial.py`, `313_anchor_allcopies.py` | partial correlations with the expression anchor |
| `242_t1_did_perm4.py` | c5 threshold sweep and the four-sample permutation |
| `310_group_contrast.py`, `311_threshold_pairwise.py` | group contrast and pairwise threshold agreement |

## Audits

| Script | Checks |
| --- | --- |
| `345_text_audit.py` | 85 numbers in the text and tables against `data/derived` |
| `346_figure_audit.py` | 91 values drawn in the figures against the text |
| `359_doi_check.py` | every reference against Crossref |
