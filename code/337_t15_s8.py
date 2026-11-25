#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
337 — 用 ROUND15 修正后的轨道重算 Table S8（11 样本两两一致，55 对）。

为什么是 11
-----------
S8 不是「r_full 有定义的样本」的集合（那是 39 个，见 S3），而是**逐拷贝面板的 11 个样本**
（4 个多能 + 7 个非多能），因此它给出的是 55 对全矩阵，而 Table 1 与 Table S3 只给
「对 H1」的一列。S8 与 Table 1 的口径必须一致，所以这里与 332 用同一批数据、同一个
SEED 与同一组染色体块重采样（共同随机数），保证同一样本对在两张表里数值相同。

Hap1 与 GM12878 本轮判为 phase_missing（未评估），涉及它们的 19 对只报 ρ_obs；
其余 36 对给 ρ_corr = ρ_obs / sqrt(r_full,i * r_full,j)，区间在同一重复内取比值。

输出：03_results/t15/t15_table_s8_pairwise.tsv
"""
import os
import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
OUT = os.path.join(RES, "t15")
HAP = os.path.join(OUT, "t15_pa5_hap_L1PA2_41.tsv.gz")
MERGED = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
UNITS = os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz")
REL = os.path.join(OUT, "t15_reliability_41.tsv")
B = 1000
SEED = 20261124
PANEL = ["H1", "H9", "WTC-11", "GM27730",
         "K562", "HepG2", "Hap1", "Hek293T", "MCF-7", "HG002", "GM12878"]


def sb_clip(r):
    if not np.isfinite(r):
        return np.nan
    return 2.0 * min(max(r, 0.0), 0.999) / (1.0 + min(max(r, 0.0), 0.999))


def main():
    hap = pd.read_csv(HAP, sep="\t")
    mrg = pd.read_csv(MERGED, sep="\t")
    units = pd.read_csv(UNITS, sep="\t")
    units = units[units.part == "L1PA_5p"].copy()
    units["grp"] = units.unit_id.str.rsplit("_", n=1).str[-1]
    units = units[units.grp == "L1PA2"].set_index("unit_id").loc[hap.unit_id].reset_index()
    assert list(units.unit_id) == list(hap.unit_id) == list(mrg.unit_id), "行序不一致"
    chrom = units.chrom.values
    rel = pd.read_csv(REL, sep="\t", keep_default_na=False, na_values=[""]).set_index("sample")

    # ---------- 块重采样（与 332 相同：SEED + 共同随机数）----------
    rng = np.random.default_rng(SEED)
    cats = pd.Categorical(chrom)
    gix = [np.where(cats.codes == j)[0] for j in range(len(cats.categories))]
    assert len(gix) == 23, len(gix)
    boots = [np.concatenate([gix[j] for j in rng.integers(0, len(gix), len(gix))])
             for _ in range(B)]

    a1 = {s: hap[f"{s}.hap1"].values for s in PANEL}
    a2 = {s: hap[f"{s}.hap2"].values for s in PANEL}
    mg = {s: mrg[s].values for s in PANEL}

    def rs(x, y, ix):
        b0 = (x[ix] == 0) & (y[ix] == 0)
        j = ix[~b0]
        if len(j) < 30 or np.ptp(x[j]) == 0 or np.ptp(y[j]) == 0:
            return np.nan
        return float(spearmanr(x[j], y[j]).statistic)

    pairs = list(combinations(PANEL, 2))
    boot_rho = {p: [] for p in pairs}
    boot_corr = {p: [] for p in pairs}
    for b, ix in enumerate(boots):
        rf = {s: sb_clip(rs(a1[s], a2[s], ix)) for s in PANEL}
        for p in pairs:
            x, y = mg[p[0]][ix], mg[p[1]][ix]
            if np.ptp(x) == 0 or np.ptp(y) == 0:
                continue
            rho = float(spearmanr(x, y).statistic)
            boot_rho[p].append(rho)
            u, v = rf[p[0]], rf[p[1]]
            if np.isfinite(u) and np.isfinite(v) and u > 0 and v > 0:
                boot_corr[p].append(rho / np.sqrt(u * v))
        if (b + 1) % 250 == 0:
            print(f"  bootstrap {b + 1}/{B}", flush=True)

    rows = []
    for i, j in pairs:
        rho = float(spearmanr(mg[i], mg[j]).statistic)
        ri = rel.loc[i, "r_full"]
        rj = rel.loc[j, "r_full"]
        ri = float(ri) if np.isfinite(ri) else np.nan
        rj = float(rj) if np.isfinite(rj) else np.nan
        corr = (rho / np.sqrt(ri * rj)) if (np.isfinite(ri) and np.isfinite(rj)
                                           and ri > 0 and rj > 0) else np.nan
        a = np.array(boot_rho[(i, j)])
        lo, hi = np.percentile(a, [2.5, 97.5]) if len(a) > 50 else (np.nan, np.nan)
        c = np.array(boot_corr[(i, j)])
        clo, chi = np.percentile(c, [2.5, 97.5]) if (len(c) > 50 and np.isfinite(corr)) else (np.nan, np.nan)
        gi = rel.loc[i, "group"]
        gj = rel.loc[j, "group"]
        gl = " x ".join(sorted([gi, gj]))
        rows.append(dict(sample_a=i, sample_b=j, group=gl, rho_obs=rho, rho_obs_lo=lo,
                         rho_obs_hi=hi, r_full_a=ri, r_full_b=rj, rho_corr=corr,
                         rho_corr_lo=clo, rho_corr_hi=chi,
                         readable_a=rel.loc[i, "readable"], readable_b=rel.loc[j, "readable"],
                         n_copy=len(mg[i])))
    s8 = pd.DataFrame(rows).sort_values(["group", "sample_a", "sample_b"])
    s8.to_csv(os.path.join(OUT, "t15_table_s8_pairwise.tsv"), sep="\t",
              index=False, na_rep="", float_format="%.4f")
    ok = s8.dropna(subset=["rho_corr"])
    print(f"\nS8: {len(s8)} 对，其中 {len(ok)} 对有 ρ_corr"
          f"（{s8.rho_corr.isna().sum()} 对涉及未评估样本）")
    print(f"ρ_corr 范围 {ok.rho_corr.min():.3f}–{ok.rho_corr.max():.3f}")
    print(ok[["sample_a", "sample_b", "rho_obs", "rho_corr"]].head(6).to_string(index=False))
    # 与 Table 1 / S3 的 H1 列逐值对齐（S8 里 H1 相关的对应等于 rho_vs_h1）
    h1 = pd.read_csv(os.path.join(OUT, "t15_rho_vs_h1_41.tsv"), sep="\t").set_index("sample")
    bad = []
    for _, r in s8.iterrows():
        if "H1" not in (r.sample_a, r.sample_b):
            continue
        o = r.sample_b if r.sample_a == "H1" else r.sample_a
        if o in h1.index and np.isfinite(h1.loc[o, "rho_obs"]):
            if abs(h1.loc[o, "rho_obs"] - r.rho_obs) > 1e-4:
                bad.append((o, h1.loc[o, "rho_obs"], r.rho_obs))
            if np.isfinite(r.rho_corr) and abs(h1.loc[o, "rho_corr"] - r.rho_corr) > 1e-4:
                bad.append((o, "corr", h1.loc[o, "rho_corr"], r.rho_corr))
    print("与 Table 1/S3 的 H1 列不一致：", bad if bad else "无")


if __name__ == "__main__":
    main()
