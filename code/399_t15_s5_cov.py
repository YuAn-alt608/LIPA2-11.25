#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
399 — Table S5（11 样本两两一致，55 对）的覆盖感知版本。

为什么重出
----------
337 的 S5 用双零口径的 r_full 做衰减校正，而 Table 1 与正文已经换成 r_full,cov：
同一对样本在两张表里给出不同的 ρ_corr（H1 0.9326 对 0.939，H1–H9 0.9833 对 0.979），
GM27730 的 0.940 也对不上 Table 1 的 0.928。S4 的同名列同理。

口径（与 389 完全一致）
----------------------
1. 逐拷贝 pa5 用覆盖碱基作分母：pa5_obs = Σv / covered_bp。
2. 一个拷贝进入 r_s 的条件是两条单倍型各覆盖 ≥800 bp。
3. r_full,cov = 2r_s,cov / (1 + r_s,cov)，不截断。
4. 区间：23 条主染色体块 bootstrap，B=1000，共同随机数（SEED 20261124），
   每个重复内重算 r_s_cov、r_full_cov 与 ρ_obs，比值在同一重复内取。
5. 只有两个样本都通过有效性规则、且 r_full,cov ≥ 0.2 时才报 ρ_corr（与 389 一致），
   并截断在 1.0。

输出：03_results/t15/t15_table_s8_pairwise_cov.tsv
"""
import os
import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import spearmanr

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")
HAP = os.path.join(T15, "t15_pa5_hap_L1PA2_41.tsv.gz")
HCOV = os.path.join(T15, "t15_hap_cov_L1PA2_41.tsv.gz")
PCOV = os.path.join(T15, "t15_pool_cov_L1PA2_41.tsv.gz")
POOL = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
UNITS = os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz")
REL = os.path.join(T15, "t15_reliability_cov_41.tsv")
OUT = os.path.join(T15, "t15_table_s8_pairwise_cov.tsv")
B = 1000
SEED = 20261124
MINBP = 800
PANEL = ["H1", "H9", "WTC-11", "GM27730",
         "K562", "HepG2", "Hap1", "Hek293T", "MCF-7", "HG002", "GM12878"]


def sb(r):
    return 2.0 * r / (1.0 + r) if np.isfinite(r) and r > -1 else np.nan


def main():
    hap = pd.read_csv(HAP, sep="\t").set_index("unit_id")
    hcov = pd.read_csv(HCOV, sep="\t").set_index("unit_id")
    pcov = pd.read_csv(PCOV, sep="\t").set_index("unit_id")
    pool = pd.read_csv(POOL, sep="\t").set_index("unit_id")
    assert list(hap.index) == list(hcov.index) == list(pcov.index) == list(pool.index)

    units = pd.read_csv(UNITS, sep="\t")
    units = units[(units.part == "L1PA_5p") & (units.unit_id.str.endswith("_L1PA2"))]
    units = units.drop_duplicates("unit_id").set_index("unit_id").loc[hap.index]
    chrom = units.chrom.values
    rel = pd.read_csv(REL, sep="\t", keep_default_na=False, na_values=[""]).set_index("sample")

    obs_h, keep, obs_pool = {}, {}, {}
    for s in PANEL:
        for h in ("hap1", "hap2"):
            p = hap[f"{s}.{h}"].values.astype(float)
            c = hcov[f"{s}.{h}"].values.astype(float)
            obs_h[(s, h)] = np.where(c > 0, p * 1000.0 / np.maximum(c, 1), np.nan)
        keep[s] = (hcov[f"{s}.hap1"].values >= MINBP) & (hcov[f"{s}.hap2"].values >= MINBP)
        p = pool[s].values.astype(float)
        c = pcov[s].values.astype(float)
        obs_pool[s] = np.where(c > 0, p * 1000.0 / np.maximum(c, 1), np.nan)

    def rs_cov(s, ix):
        k = keep[s][ix] & np.isfinite(obs_h[(s, "hap1")][ix]) & np.isfinite(obs_h[(s, "hap2")][ix])
        j = ix[k]
        if len(j) < 30:
            return np.nan
        a, b = obs_h[(s, "hap1")][j], obs_h[(s, "hap2")][j]
        if np.ptp(a) == 0 or np.ptp(b) == 0:
            return np.nan
        return float(spearmanr(a, b).statistic)

    # ---------- 块重采样（与 389/337 同一构造）----------
    cats = pd.Categorical(chrom)
    gix = [np.where(cats.codes == j)[0] for j in range(len(cats.categories))]
    assert len(gix) == 23, len(gix)
    rng = np.random.default_rng(SEED)
    boots = [np.concatenate([gix[j] for j in rng.integers(0, len(gix), len(gix))])
             for _ in range(B)]

    pairs = list(combinations(PANEL, 2))
    boot_rho = {p: [] for p in pairs}
    boot_corr = {p: [] for p in pairs}
    for b, ix in enumerate(boots):
        rf = {s: sb(rs_cov(s, ix)) for s in PANEL}
        for p in pairs:
            x, y = obs_pool[p[0]][ix], obs_pool[p[1]][ix]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() <= 30 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
                continue
            rho = float(spearmanr(x[ok], y[ok]).statistic)
            boot_rho[p].append(rho)
            u, v = rf[p[0]], rf[p[1]]
            if np.isfinite(u) and np.isfinite(v) and u > 0 and v > 0:
                boot_corr[p].append(min(rho / np.sqrt(u * v), 1.0))
        if (b + 1) % 250 == 0:
            print(f"  bootstrap {b + 1}/{B}", flush=True)

    rows = []
    for i, j in pairs:
        x, y = obs_pool[i], obs_pool[j]
        ok = np.isfinite(x) & np.isfinite(y)
        rho = float(spearmanr(x[ok], y[ok]).statistic)
        ri, rj = rel.loc[i, "r_full_cov"], rel.loc[j, "r_full_cov"]
        ri = float(ri) if np.isfinite(ri) else np.nan
        rj = float(rj) if np.isfinite(rj) else np.nan
        vi, vj = bool(rel.loc[i, "valid"]), bool(rel.loc[j, "valid"])
        corr = (min(rho / np.sqrt(ri * rj), 1.0)
                if (np.isfinite(ri) and np.isfinite(rj) and ri >= 0.2 and rj >= 0.2 and vi and vj)
                else np.nan)
        a = np.array(boot_rho[(i, j)])
        lo, hi = np.percentile(a, [2.5, 97.5]) if len(a) > 50 else (np.nan, np.nan)
        c = np.array(boot_corr[(i, j)])
        clo, chi = np.percentile(c, [2.5, 97.5]) if (len(c) > 50 and np.isfinite(corr)) else (np.nan, np.nan)
        gl = " x ".join(sorted([rel.loc[i, "group"], rel.loc[j, "group"]]))
        rows.append(dict(sample_a=i, sample_b=j, group=gl, rho_obs=rho, rho_obs_lo=lo,
                         rho_obs_hi=hi,
                         r_s_cov_a=rel.loc[i, "r_s_cov"], r_full_cov_a=ri,
                         r_s_cov_b=rel.loc[j, "r_s_cov"], r_full_cov_b=rj,
                         rho_corr=corr, rho_corr_lo=clo, rho_corr_hi=chi,
                         readable_a=rel.loc[i, "readable"], readable_b=rel.loc[j, "readable"],
                         n_copy=int(ok.sum())))
    s5 = pd.DataFrame(rows).sort_values(["group", "sample_a", "sample_b"])
    s5.to_csv(OUT, sep="\t", index=False, na_rep="", float_format="%.4f")

    ok = s5.dropna(subset=["rho_corr"])
    print(f"\nS5(cov): {len(s5)} 对，其中 {len(ok)} 对有 ρ_corr")
    print(f"ρ_corr 范围 {ok.rho_corr.min():.3f}–{ok.rho_corr.max():.3f}")
    pl = ok[ok.group.str.contains("pluripotent x pluripotent")]
    print(f"多能对 {len(pl)}：{pl.rho_corr.min():.3f}–{pl.rho_corr.max():.3f}")
    npl = ok[ok.group.str.contains("pluripotent") & ~ok.group.str.contains("pluripotent x pluripotent")]
    print(f"多能-非多能对 {len(npl)}：{npl.rho_corr.min():.3f}–{npl.rho_corr.max():.3f}")
    print(f"H1 对：{ok[(ok.sample_a == 'H1') | (ok.sample_b == 'H1')].rho_corr.min():.3f}"
          f"–{ok[(ok.sample_a == 'H1') | (ok.sample_b == 'H1')].rho_corr.max():.3f}")

    # 与 Table 1（对 H1 的一列）逐值对齐
    h1 = pd.read_csv(os.path.join(T15, "t15_rho_vs_h1_cov_41.tsv"), sep="\t").set_index("sample")
    bad = []
    for _, r in s5.iterrows():
        if "H1" not in (r.sample_a, r.sample_b):
            continue
        o = r.sample_b if r.sample_a == "H1" else r.sample_a
        if o in h1.index and np.isfinite(h1.loc[o, "rho_obs"]):
            if abs(h1.loc[o, "rho_obs"] - r.rho_obs) > 1e-3:
                bad.append((o, "obs", h1.loc[o, "rho_obs"], round(r.rho_obs, 4)))
            if np.isfinite(r.rho_corr) and abs(h1.loc[o, "rho_corr"] - r.rho_corr) > 1e-3:
                bad.append((o, "corr", h1.loc[o, "rho_corr"], round(r.rho_corr, 4)))
    print("与 Table 1 的 H1 列不一致：", bad if bad else "无")


if __name__ == "__main__":
    main()
