#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
333 — ROUND15 的增补产物：阈值敏感性 S7、定相诊断（Fig. S3）、Fig. 2 面板数据、
以及「多能样本 ρ_corr 0.94–0.98」这一句的复核。

全部从 t15_reliability_41.tsv 与修正后的 hap 表读数，不引用任何手抄值。

产出 03_results/t15/
  t15_table_s7_threshold.tsv   阈值 0.5–0.9，分母 40（通过定相检查的样本）
  t15_selfconsistency.tsv      每样本 (hap1+hap2)/2 vs 自身 pooled —— 抓出本轮回题的量
  t15_phasing_summary.tsv      计数汇总（41 / 40 / 10 / 30 / 1）
  t15_fig2_panel_data.tsv      Fig. 2a/2b 的逐样本数据
  t15_pluri_pairwise.tsv       4 个多能样本两两 ρ_obs / ρ_corr（复核 0.94–0.98）
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
OUT = os.path.join(RES, "t15")
HAP = os.path.join(OUT, "t15_pa5_hap_L1PA2_41.tsv.gz")
MERGED = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
REL = os.path.join(OUT, "t15_reliability_41.tsv")
THR = [0.5, 0.6, 0.667, 0.7, 0.75, 0.8, 0.9]
PLURI = ["H1", "H9", "WTC-11", "GM27730"]


def rho(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10 or np.ptp(x[m]) == 0 or np.ptp(y[m]) == 0:
        return np.nan
    return float(spearmanr(x[m], y[m]).statistic)


def main():
    # readable 含字面量 "n/a"，pandas 默认会当缺失值读成 NaN，必须关掉
    rel = pd.read_csv(REL, sep="\t", keep_default_na=False,
                      na_values=[""])
    hap = pd.read_csv(HAP, sep="\t")
    mrg = pd.read_csv(MERGED, sep="\t")
    assert len(rel) == 41

    ph = rel[rel["phased"].astype(bool)].copy()
    rd = rel[rel["readable"].isin(["yes", "borderline"])]

    # ---- S7 ----
    rows = []
    for t in THR:
        above = ph[ph["r_full_lo"] >= t]
        span = ph[(ph["r_full_lo"] < t) & (ph["r_full_hi"] >= t)]
        rows.append(dict(
            threshold=t, n_phased=len(ph),
            readable=int((ph["r_full"] >= t).sum()),
            ci_entirely_above=len(above), ci_spans_threshold=len(span),
            below_point=int((ph["r_full"] < t).sum()),
            ci_entirely_below=int((ph["r_full_hi"] < t).sum()),
            readable_samples=", ".join(sorted(ph.loc[ph["r_full"] >= t, "sample"])),
            borderline_samples=", ".join(sorted(span["sample"])),
        ))
    s7 = pd.DataFrame(rows)
    s7.to_csv(os.path.join(OUT, "t15_table_s7_threshold.tsv"), sep="\t", index=False)
    print("S7\n", s7.drop(columns=["readable_samples", "borderline_samples"]).to_string(index=False))

    # ---- 自洽性：(hap1+hap2)/2 vs 自身 pooled ----
    rows = []
    for s in rel["sample"]:
        a1, a2, m = hap[f"{s}.hap1"].values, hap[f"{s}.hap2"].values, mrg[s].values
        rows.append(dict(sample=s, group=rel.loc[rel["sample"] == s, "group"].iloc[0],
                         phased=bool(rel.loc[rel["sample"] == s, "phased"].iloc[0]),
                         rho_mean_vs_pooled=round(rho((a1 + a2) / 2, m), 4),
                         rho_h1_vs_pooled=round(rho(a1, m), 4),
                         rho_h2_vs_pooled=round(rho(a2, m), 4)))
    sc = pd.DataFrame(rows)
    sc.to_csv(os.path.join(OUT, "t15_selfconsistency.tsv"), sep="\t", index=False)
    g = sc[sc.phased]
    print(f"\n自洽性：{len(g)} 个通过定相的样本 (hap1+hap2)/2 vs 自身 pooled "
          f"{g.rho_mean_vs_pooled.min():.3f}–{g.rho_mean_vs_pooled.max():.3f}"
          f"（中位 {g.rho_mean_vs_pooled.median():.3f}）")

    # ---- 汇总 ----
    summ = pd.DataFrame([
        dict(metric="samples_analysed", value=len(rel)),
        dict(metric="not_assessed_frac_both0", value=int((rel["note"].fillna("")
                                                         .str.contains("frac_both0")).sum())),
        dict(metric="not_assessed_rho", value=int((rel["note"].fillna("")
                                                  .str.contains("rho\\(merged")).sum())),
        dict(metric="phased", value=len(ph)),
        dict(metric="readable_yes", value=int((rel["readable"] == "yes").sum())),
        dict(metric="readable_borderline", value=int((rel["readable"] == "borderline").sum())),
        dict(metric="not_readable", value=int((rel["readable"] == "no").sum())),
        dict(metric="r_s_nz_le0_but_phased", value=int(rel["r_s_le0"].sum())),
        dict(metric="rho_pooled_hap_min_phased",
             value=round(float(np.nanmin(ph[["rho_m_h1", "rho_m_h2"]].min(axis=1))), 4)),
        dict(metric="r_full_max_not_readable",
             value=round(float(rel.loc[rel["readable"] == "no", "r_full"].max()), 4)),
    ])
    summ.to_csv(os.path.join(OUT, "t15_phasing_summary.tsv"), sep="\t", index=False)
    print("\n", summ.to_string(index=False))

    # ---- Fig. 2a/2b ----
    fig = rel[["sample", "group", "med_pa5", "r_s_all", "r_s_nz", "r_full", "r_full_lo",
               "r_full_hi", "readable", "phased", "frac_both0"]].copy()
    fig.to_csv(os.path.join(OUT, "t15_fig2_panel_data.tsv"), sep="\t", index=False)

    # ---- 多能样本两两 ----
    rows = []
    for i, a in enumerate(PLURI):
        for b in PLURI[i + 1:]:
            r = rho(mrg[a].values, mrg[b].values)
            fa = float(rel.loc[rel["sample"] == a, "r_full"].iloc[0])
            fb = float(rel.loc[rel["sample"] == b, "r_full"].iloc[0])
            rows.append(dict(sample_a=a, sample_b=b, rho_obs=round(r, 4),
                             r_full_a=fa, r_full_b=fb,
                             rho_corr=round(r / np.sqrt(fa * fb), 4)))
    pw = pd.DataFrame(rows)
    pw.to_csv(os.path.join(OUT, "t15_pluri_pairwise.tsv"), sep="\t", index=False)
    print(f"\n多能两两 ρ_corr：{pw.rho_corr.min():.3f}–{pw.rho_corr.max():.3f}")
    print(pw.to_string(index=False))

    # ---- 可读名单 ----
    print("\n可读:", ", ".join(sorted(rd["sample"])))
    print("不可读:", ", ".join(sorted(rel.loc[rel["readable"] == "no", "sample"])))
    print("未评估:", ", ".join(sorted(rel.loc[rel["readable"] == "n/a", "sample"])))


if __name__ == "__main__":
    main()
