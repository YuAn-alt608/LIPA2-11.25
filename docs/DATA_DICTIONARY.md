# Data dictionary

Every file below is a frozen copy of the analysis table used by the manuscript, 
the figures or Tables S1–S8. Paths are relative to `data/derived/`. 
All tables are tab-separated; `.gz` files are gzip-compressed TSV. 
Row counts exclude the header.

## `l1pa5_copy_tracks_v2.tsv.gz`

Rows: 14444. Columns: 12.

| column |
| --- |
| `unit_id` |
| `part` |
| `chrom` |
| `start` |
| `end` |
| `h1_k9` |
| `h1_dnase` |
| `h9_k9` |
| `h9_acc` |
| `gc` |
| `cpg_oe` |
| `h1_pm` |

## `sample_labels.tsv`

Rows: 41. Columns: 7.

| column |
| --- |
| `sample` |
| `cell_class` |
| `pluripotent` |
| `source` |
| `cellosaurus` |
| `cs_category` |
| `cs_verdict` |

## `t10/t10_readability_map.tsv`

Rows: 143. Columns: 19.

| column |
| --- |
| `family` |
| `sample` |
| `pa5_median` |
| `r_s` |
| `note` |
| `r_full` |
| `n_hap` |
| `cov_med` |
| `frac_nonzero` |
| `source` |
| `r_s_all` |
| `r_s_nz` |
| `r_full_all` |
| `frac_both0` |
| `rho_merged_hap1` |
| `rho_merged_hap2` |
| `med_pa5_both0` |
| `med_pa5_rest` |
| `phase_valid` |

## `t10_sim_band.tsv`

Rows: 13. Columns: 15.

| column |
| --- |
| `sample` |
| `obs_rho` |
| `band_lo` |
| `band_hi` |
| `sim_med` |
| `in_band` |
| `direction` |
| `sigma` |
| `r_full` |
| `band_lo_2side` |
| `band_hi_2side` |
| `in_band_2side` |
| `sigma_H1` |
| `pa5_median` |
| `is_pluri` |

## `t13/t13_group_contrast_sample_first.tsv`

Rows: 30. Columns: 15.

| column |
| --- |
| `family` |
| `window` |
| `n_copy` |
| `group_pluri` |
| `group_other` |
| `raw_diff` |
| `raw_lo` |
| `raw_hi` |
| `ref_raw_diff` |
| `ref_lo` |
| `ref_hi` |
| `did` |
| `did_lo` |
| `did_hi` |
| `reference` |

## `t13/t13_partial_rho_allcopies.tsv`

Rows: 16. Columns: 9.

| column |
| --- |
| `scope` |
| `pred` |
| `n` |
| `n_cov` |
| `n_amb` |
| `rho_partial` |
| `lo` |
| `hi` |
| `rho_raw` |

## `t13/t13_per_sample_did.tsv`

Rows: 11. Columns: 7.

| column |
| --- |
| `sample` |
| `group` |
| `l1pa2_median` |
| `l1pa4_7_median` |
| `diff` |
| `did_after_loo` |
| `role` |

## `t15/t15_pa5_hap_L1PA2_41.tsv.gz`

Rows: 1032. Columns: 83.

| column |
| --- |
| `unit_id` |
| `GM12878.hap1` |
| `GM12878.hap2` |
| `K562.hap1` |
| `K562.hap2` |
| `HepG2.hap1` |
| `HepG2.hap2` |
| `H1.hap1` |
| `H1.hap2` |
| `Hap1.hap1` |
| `Hap1.hap2` |
| `Hek293T.hap1` |
| `Hek293T.hap2` |
| `Jurkat.hap1` |
| `Jurkat.hap2` |
| `HG002.hap1` |
| `HG002.hap2` |
| `Panc1.hap1` |
| `Panc1.hap2` |
| `THP-1.hap1` |
| `THP-1.hap2` |
| `H9.hap1` |
| `H9.hap2` |
| `HG01123.hap1` |
| `HG01123.hap2` |
| `HG01258.hap1` |
| `HG01258.hap2` |
| `HG01358.hap1` |
| `HG01358.hap2` |
| `HG01361.hap1` |
| `HG01361.hap2` |
| `HG01891.hap1` |
| `HG01891.hap2` |
| `HG02071.hap1` |
| `HG02071.hap2` |
| `HG02074.hap1` |
| `HG02074.hap2` |
| `HG02132.hap1` |
| `HG02132.hap2` |
| `HG02135.hap1` |
| `HG02135.hap2` |
| `HG02257.hap1` |
| `HG02257.hap2` |
| `HG02486.hap1` |
| `HG02486.hap2` |
| `HG02559.hap1` |
| `HG02559.hap2` |
| `HG02572.hap1` |
| `HG02572.hap2` |
| `HG02717.hap1` |
| `HG02717.hap2` |
| `HG02886.hap1` |
| `HG02886.hap2` |
| `HG03516.hap1` |
| `HG03516.hap2` |
| `HG03804.hap1` |
| `HG03804.hap2` |
| `HG03942.hap1` |
| `HG03942.hap2` |
| `HG04160.hap1` |
| `HG04160.hap2` |
| `HG04187.hap1` |
| `HG04187.hap2` |
| `A549.hap1` |
| `A549.hap2` |
| `MCF-7.hap1` |
| `MCF-7.hap2` |
| `HCT116.hap1` |
| `HCT116.hap2` |
| `Caco-2.hap1` |
| `Caco-2.hap2` |
| `WTC-11.hap1` |
| `WTC-11.hap2` |
| `GM28572.hap1` |
| `GM28572.hap2` |
| `GM28570.hap1` |
| `GM28570.hap2` |
| `GM25456.hap1` |
| `GM25456.hap2` |
| `GM25455.hap1` |
| `GM25455.hap2` |
| `GM27730.hap1` |
| `GM27730.hap2` |

## `t15/t15_pa5_vs_rfull.tsv`

Rows: 6. Columns: 4.

| column |
| --- |
| `subset` |
| `n` |
| `rho` |
| `p` |

## `t15/t15_reliability_41.tsv`

Rows: 41. Columns: 19.

| column |
| --- |
| `sample` |
| `group` |
| `n` |
| `n_nz` |
| `frac_both0` |
| `nz_hap1` |
| `nz_hap2` |
| `med_pa5` |
| `r_s_all` |
| `r_s_nz` |
| `rho_m_h1` |
| `rho_m_h2` |
| `phased` |
| `note` |
| `r_full` |
| `r_full_lo` |
| `r_full_hi` |
| `readable` |
| `r_s_le0` |

## `t15/t15_rho_vs_h1_41.tsv`

Rows: 40. Columns: 8.

| column |
| --- |
| `sample` |
| `rho_obs` |
| `r_full` |
| `rho_corr` |
| `lo` |
| `hi` |
| `n_boot` |
| `readable` |

## `t15/t15_selfconsistency.tsv`

Rows: 41. Columns: 6.

| column |
| --- |
| `sample` |
| `group` |
| `phased` |
| `rho_mean_vs_pooled` |
| `rho_h1_vs_pooled` |
| `rho_h2_vs_pooled` |

## `t15/t15_table_s3.tsv`

Rows: 41. Columns: 23.

| column |
| --- |
| `sample` |
| `group` |
| `batch` |
| `n` |
| `n_nz` |
| `frac_both0` |
| `nz_hap1` |
| `nz_hap2` |
| `med_pa5` |
| `r_s_all` |
| `r_s_nz` |
| `rho_m_h1` |
| `rho_m_h2` |
| `phased` |
| `note` |
| `r_full` |
| `r_full_lo` |
| `r_full_hi` |
| `readable` |
| `rho_obs` |
| `rho_corr` |
| `lo` |
| `hi` |

## `t15/t15_table_s7_threshold.tsv`

Rows: 7. Columns: 9.

| column |
| --- |
| `threshold` |
| `n_phased` |
| `readable` |
| `ci_entirely_above` |
| `ci_spans_threshold` |
| `below_point` |
| `ci_entirely_below` |
| `readable_samples` |
| `borderline_samples` |

## `t15/t15_table_s8_pairwise.tsv`

Rows: 55. Columns: 14.

| column |
| --- |
| `sample_a` |
| `sample_b` |
| `group` |
| `rho_obs` |
| `rho_obs_lo` |
| `rho_obs_hi` |
| `r_full_a` |
| `r_full_b` |
| `rho_corr` |
| `rho_corr_lo` |
| `rho_corr_hi` |
| `readable_a` |
| `readable_b` |
| `n_copy` |

## `t1_did_by_sample.tsv`

Rows: 41. Columns: 12.

| column |
| --- |
| `sample` |
| `open5_L1PA2` |
| `open5_old_mean` |
| `did` |
| `med5_L1PA2` |
| `pluripotent` |
| `n_peaks` |
| `bp_covered` |
| `genome_fire_fraction` |
| `log_depth` |
| `measure` |
| `cell_class` |

## `t1_did_perm4_c5.tsv`

Rows: 11. Columns: 21.

| column |
| --- |
| `threshold` |
| `top4` |
| `top4_is_group` |
| `group_mean` |
| `best_other` |
| `ratio` |
| `perm_p` |
| `n_comb` |
| `rank_H1` |
| `rank_H9` |
| `rank_WTC11` |
| `rank_GM27730` |
| `intruder` |
| `old_top3` |
| `old_top3_is_group` |
| `old_perm_p` |
| `old_n_comb` |
| `pluri_open_copies_median` |
| `in_prespecified_range` |
| `measure` |
| `measure_note` |

## `t1_did_perm4_pa5.tsv`

Rows: 1. Columns: 23.

| column |
| --- |
| `subset` |
| `n_units` |
| `top4` |
| `top4_is_group` |
| `group_mean` |
| `best_other` |
| `best_other_name` |
| `ratio` |
| `perm_p` |
| `n_comb` |
| `rank_GM27730` |
| `top4_rank4_value` |
| `H1` |
| `H9` |
| `WTC-11` |
| `K562` |
| `HepG2` |
| `GM27730` |
| `Hap1` |
| `Hek293T` |
| `MCF-7` |
| `measure` |
| `measure_note` |

## `t1d_pa5_by_copy_9.tsv.gz`

Rows: 7222. Columns: 13.

| column |
| --- |
| `unit_id` |
| `grp` |
| `H1` |
| `H9` |
| `WTC-11` |
| `K562` |
| `HepG2` |
| `GM27730` |
| `Hap1` |
| `Hek293T` |
| `MCF-7` |
| `measure` |
| `measure_note` |

## `t9b/t9b_pa5_L1PA2_41_fixed.tsv.gz`

Rows: 1032. Columns: 42.

| column |
| --- |
| `unit_id` |
| `GM12878` |
| `K562` |
| `HepG2` |
| `H1` |
| `Hap1` |
| `Hek293T` |
| `Jurkat` |
| `HG002` |
| `Panc1` |
| `THP-1` |
| `H9` |
| `HG01123` |
| `HG01258` |
| `HG01358` |
| `HG01361` |
| `HG01891` |
| `HG02071` |
| `HG02074` |
| `HG02132` |
| `HG02135` |
| `HG02257` |
| `HG02486` |
| `HG02559` |
| `HG02572` |
| `HG02717` |
| `HG02886` |
| `HG03516` |
| `HG03804` |
| `HG03942` |
| `HG04160` |
| `HG04187` |
| `A549` |
| `MCF-7` |
| `HCT116` |
| `Caco-2` |
| `WTC-11` |
| `GM28572` |
| `GM28570` |
| `GM25456` |
| `GM25455` |
| `GM27730` |

## `t9d/t9d_coupling.tsv`

Rows: 6. Columns: 7.

| column |
| --- |
| `meth_track` |
| `pa5_of` |
| `rho` |
| `lo` |
| `hi` |
| `n` |
| `same_cells` |

## `t9i/t10_dA_dM_controls.tsv`

Rows: 9. Columns: 8.

| column |
| --- |
| `control` |
| `stratum` |
| `n` |
| `rho` |
| `lo` |
| `hi` |
| `p` |
| `note` |

## `t9i/t10_dA_dM_families.tsv`

Rows: 3. Columns: 5.

| column |
| --- |
| `family` |
| `n` |
| `rho` |
| `lo` |
| `hi` |

## `t9i/t10_dA_dM_main.tsv`

Rows: 1. Columns: 24.

| column |
| --- |
| `n` |
| `rho` |
| `lo` |
| `hi` |
| `p_perm` |
| `null_lo` |
| `null_hi` |
| `concord` |
| `p_binom` |
| `E1` |
| `E2` |
| `E3` |
| `sd_dA` |
| `sd_dM` |
| `med_dA` |
| `med_dM` |
| `rho_common` |
| `rho_partial` |
| `r_dcp_dM` |
| `r_dcp_dA` |
| `frac_cpg_set_identical` |
| `med_only_one_cpg` |
| `n_cpg1_med` |
| `n_cpg2_med` |

## `t9i/t10_dA_dM_quintile.tsv`

Rows: 5. Columns: 4.

| column |
| --- |
| `quintile` |
| `n` |
| `rho` |
| `abs_dA_med` |

## `t9i/t10_dA_dM_regions.tsv`

Rows: 4. Columns: 5.

| column |
| --- |
| `region` |
| `n` |
| `rho` |
| `lo` |
| `hi` |

## `t9i/t9i_dA_by_copy.tsv.gz`

Rows: 1032. Columns: 14.

| column |
| --- |
| `unit_id` |
| `chrom` |
| `start` |
| `end` |
| `pa5_GM27730` |
| `pa5_HG002` |
| `dA` |
| `pa5_H1` |
| `dA_hap_H1` |
| `dM_hap_H1` |
| `dM_common_cpg` |
| `nCpg_hap1` |
| `nCpg_hap2` |
| `nCpg_union` |

## `t9m/t9m_partial_rho.tsv`

Rows: 16. Columns: 8.

| column |
| --- |
| `scope` |
| `pred` |
| `n` |
| `n_cov` |
| `rho_partial` |
| `lo` |
| `hi` |
| `rho_raw` |

## `t9m/t9m_rna_by_copy.tsv.gz`

Rows: 1032. Columns: 36.

| column |
| --- |
| `unit_id` |
| `grp` |
| `chrom` |
| `el_start` |
| `el_end` |
| `start` |
| `end` |
| `strand` |
| `len` |
| `cat` |
| `cgi` |
| `tss_dist` |
| `gc` |
| `cpg_oe` |
| `map50` |
| `label_body` |
| `label_5p` |
| `senseb_r1` |
| `senseb_r2` |
| `senseb_r3` |
| `antib_r1` |
| `antib_r2` |
| `antib_r3` |
| `sense5_r1` |
| `sense5_r2` |
| `sense5_r3` |
| `anti5_r1` |
| `anti5_r2` |
| `anti5_r3` |
| `pa5_H1` |
| `pa5_H9` |
| `pa5_WTC-11` |
| `pa5_GM27730` |
| `pa5_K562` |
| `pa5_HepG2` |
| `pa5_HG002` |
