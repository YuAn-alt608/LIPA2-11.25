#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
391 — Table S4 的 r28 版本：覆盖碱基分母 + 80% 覆盖门槛。

与 375（r26/r27）的三处差别
---------------------------
1. 观测值改用分母修正后的 pa5_obs = Σv / covered_bp（原为 Σv / 1,000）。
2. 保留规则从「两 hap 各覆盖 ≥1 bp」改为「各 ≥800 bp」。r_s_cov / r_full_cov 随
   之重算；门槛从 1 bp 提到 800 bp 的敏感性单独列出（r_s_cov_1bp）。
3. 新增 med_frac_cov：该样本池化轨道在这 1,032 个拷贝上的覆盖率中位数，用于把
   「27.8% 是合并数字、2–92% 是单样本范围」拆开（评审要点）。

输出（覆盖 375 的同名文件，先备份）
----
t15/t15_table_s4_cov.tsv
t15/t15_table_s4_cov_summary.tsv
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")
HAP = os.path.join(T15, "t15_pa5_hap_L1PA2_41.tsv.gz")
HCOV = os.path.join(T15, "t15_hap_cov_L1PA2_41.tsv.gz")
PCOV = os.path.join(T15, "t15_pool_cov_L1PA2_41.tsv.gz")
POOL = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
S3 = os.path.join(T15, "t15_table_s3.tsv")
OUT = os.path.join(T15, "t15_table_s4_cov.tsv")
SUM = os.path.join(T15, "t15_table_s4_cov_summary.tsv")
MINBP = 800
B = 1000
THR = 0.667


def sb(r):
    return max(0.0, 2 * r / (1 + r)) if np.isfinite(r) and r > -1 else np.nan


def main():
    hap = pd.read_csv(HAP, sep="\t").set_index("unit_id")
    cv = pd.read_csv(HCOV, sep="\t").set_index("unit_id")
    pc = pd.read_csv(PCOV, sep="\t").set_index("unit_id")
    pool = pd.read_csv(POOL, sep="\t").set_index("unit_id")
    s4 = pd.read_csv(S3, sep="\t")
    samples = [c[:-5] for c in hap.columns if c.endswith(".hap1")]

    rows = {}
    tot_both0 = tot_cov = 0
    for s in samples:
        a, b = hap[f"{s}.hap1"].values, hap[f"{s}.hap2"].values
        ca, cb = cv[f"{s}.hap1"].values, cv[f"{s}.hap2"].values
        oa = np.where(ca > 0, a * B / np.maximum(ca, 1), np.nan)
        ob = np.where(cb > 0, b * B / np.maximum(cb, 1), np.nan)
        both0 = (a == 0) & (b == 0)
        keep80 = np.isfinite(oa) & np.isfinite(ob) & (ca >= MINBP) & (cb >= MINBP)
        keep1 = (ca > 0) & (cb > 0)
        n_both0 = int(both0.sum())
        n_both0_cov = int((both0 & (ca >= MINBP) & (cb >= MINBP)).sum())
        r80 = spearmanr(oa[keep80], ob[keep80]).statistic if keep80.sum() > 5 else np.nan
        r1 = spearmanr(oa[keep1], ob[keep1]).statistic if keep1.sum() > 5 else np.nan
        pfrac = pc[s].values / B
        rows[s] = dict(n_cov=int(keep80.sum()), frac_both0_covered=round(n_both0_cov / n_both0, 4) if n_both0 else np.nan,
                       r_s_cov=None if not np.isfinite(r80) else round(float(r80), 4),
                       r_full_cov=None if not np.isfinite(sb(r80)) else round(float(sb(r80)), 4),
                       r_s_cov_1bp=None if not np.isfinite(r1) else round(float(r1), 4),
                       med_frac_cov=round(float(np.nanmedian(pfrac)), 3),
                       n_both0=n_both0, n_both0_cov800=n_both0_cov)
        tot_both0 += n_both0
        tot_cov += n_both0_cov

    d = s4.copy()
    for k in ("n_cov", "frac_both0_covered", "r_s_cov", "r_full_cov", "r_s_cov_1bp",
              "med_frac_cov", "n_both0_cov800"):
        if k in d.columns:
            d = d.drop(columns=[k])
    add = pd.DataFrame(rows).T.rename_axis("sample").reset_index()
    d = d.merge(add, on="sample", how="left")
    # readable 改为按 r_full,cov 判定（旧列来自 r27 的双零口径）
    def _read(row):
        if not bool(row["phased"]) or not np.isfinite(row["r_full_cov"] or np.nan):
            return "n/a"
        if row["r_full_cov"] < THR:
            return "no"
        return "borderline" if row["r_full_cov_lo"] < THR else "yes"
    lo = pd.read_csv("../03_results/t15/t15_reliability_cov_41.tsv", sep="\t").set_index("sample")
    d = d.merge(lo[["r_full_cov_lo", "r_full_cov_hi", "med_pa5_cov", "frac_nocov", "valid"]],
                on="sample", how="left")
    d["readable"] = d.apply(_read, axis=1)
    # phased / note 按正文的有效性规则重算（两条单倍型都没有覆盖），不再沿用 S3 的双零判定
    d["phased"] = d["valid"].astype(bool)
    d["note"] = np.where(d["frac_nocov"] > 0.5, "both haplotypes uncovered in >50% of copies",
                         np.where(d[["rho_m_h1", "rho_m_h2"]].min(axis=1) < 0.2,
                                  "merged-hap rho<0.2", ""))
    d = d.drop(columns=["valid"])
    d.to_csv(OUT, sep="\t", index=False)

    a = d[d.n_cov.notna() & d.phased]
    nr = a[a.readable == "no"]
    rd = a[a.readable.isin(["yes", "borderline"])]
    lcl = a[a.group == "other (LCL)"]
    summ = pd.DataFrame([dict(
        n_with_track=int(d.n_cov.notna().sum()), n_assessable=len(a),
        tot_both0=int(tot_both0), tot_both0_cov=int(tot_cov),
        pct_both0_covered=round(100 * tot_cov / tot_both0, 2),
        frac_min=round(float(a.frac_both0_covered.min()), 3),
        frac_max=round(float(a.frac_both0_covered.max()), 3),
        frac_median=round(float(a.frac_both0_covered.median()), 3),
        medfrac_min=round(float(a.med_frac_cov.min()), 3),
        medfrac_median=round(float(a.med_frac_cov.median()), 3),
        medfrac_max=round(float(a.med_frac_cov.max()), 3),
        r_s_min=round(float(a.r_s_cov.min()), 4), r_s_max=round(float(a.r_s_cov.max()), 4),
        n_le0=int((a.r_s_cov <= 0).sum()), n_ge_half=int((a.r_s_cov >= 0.5).sum()),
        r_full_max=round(float(a.r_full_cov.max()), 3),
        nr_r_s_min=round(float(nr.r_s_cov.min()), 4), nr_r_s_max=round(float(nr.r_s_cov.max()), 4),
        rd_r_s_min=round(float(rd.r_s_cov.min()), 4), rd_r_s_max=round(float(rd.r_s_cov.max()), 4),
        lcl_r_s_min=round(float(lcl.r_s_cov.min()), 4), lcl_r_s_max=round(float(lcl.r_s_cov.max()), 4),
        n_readable=int(a.readable.isin(["yes", "borderline"]).sum()))])
    summ.to_csv(SUM, sep="\t", index=False)
    print(summ.T.to_string())
    print()
    print(d[["sample", "n_cov", "med_frac_cov", "frac_both0_covered", "r_s_cov",
             "r_full_cov", "r_s_nz", "r_full", "readable"]].to_string(index=False))


if __name__ == "__main__":
    main()
