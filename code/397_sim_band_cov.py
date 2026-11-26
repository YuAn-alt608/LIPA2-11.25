#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
397 — 噪声带改用 r_full,cov 重算。

303 的 t10_sim_band.tsv 用的是双零口径的 r_full（11 月 25 日）。图 3b、图 S4 与
"H9 高出噪声带 +0.024" 都读它，而 ρ_corr 已经换成 r_full,cov，输入不一致。

这里只重算噪声带，公式与 303 完全一致：
  σ      = sqrt(var(mapped) * (1/r_full - 1))         目标样本的噪声
  σ_H1   = sqrt(var(pa5_H1) * (1/r_full_H1 - 1))      H1 的噪声
  band   = H1 排序映射到目标样本分布后加噪，再与 H1 求 ρ（单侧噪声）
  band2  = 两侧都加噪（敏感性）
r_full 取 t15_reliability_cov_41.tsv 的 r_full_cov，H1 取同表的 0.939。

输出 t10/t10_sim_band_cov.tsv，列与 t10_sim_band.tsv 兼容。
"""
import os
import numpy as np
import pandas as pd

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T9B = os.path.join(RES, "t9b")
T15 = os.path.join(RES, "t15")
T10 = os.path.join(RES, "t10")
PA_FIXED = os.path.join(T9B, "t9b_pa5_L1PA2_41_fixed.tsv.gz")
REL = os.path.join(T15, "t15_reliability_cov_41.tsv")
FULL = os.path.join(T15, "t15_reliability_41.tsv")
S4 = os.path.join(T15, "t15_table_s4_cov.tsv")
OUT = os.path.join(RES, "t10_sim_band_cov.tsv")
B_SIM = 200
SEED = 20261124
PLURI = ["H1", "H9", "WTC-11", "GM27730"]


def sp(x, y):
    x = pd.Series(np.asarray(x, float))
    y = pd.Series(np.asarray(y, float))
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan
    return float(x[m].rank().corr(y[m].rank()))


def main():
    pa = pd.read_csv(PA_FIXED, sep="\t")
    rel = pd.read_csv(REL, sep="\t").set_index("sample")
    old = pd.read_csv(FULL, sep="\t").set_index("sample")
    s4 = pd.read_csv(S4, sep="\t").set_index("sample")
    pa = pa.set_index("unit_id")
    paH1 = pa["H1"].to_numpy(float)
    rfull_H1 = float(rel.loc["H1", "r_full_cov"])

    rows = []
    for s in pa.columns:
        if s == "H1":
            continue
        if s not in rel.index or not np.isfinite(rel.loc[s, "r_full_cov"]):
            continue
        r_full = float(rel.loc[s, "r_full_cov"])
        tgt = pa[s].to_numpy(float)
        h9 = pa["H9"].to_numpy(float)
        m = np.isfinite(h9) & np.isfinite(tgt)
        if m.sum() < 30:
            continue
        qtl = (pd.Series(h9[m]).rank().to_numpy() - 0.5) / m.sum()
        mapped = np.quantile(tgt[m], qtl)
        var = float(np.var(mapped)) or 1.0
        sigma = np.sqrt(max(var * (1.0 / max(r_full, 0.05) - 1.0), 1e-9))
        varH1 = float(np.var(paH1[m])) or 1.0
        sigmaH1 = np.sqrt(max(varH1 * (1.0 / max(rfull_H1, 0.05) - 1.0), 1e-9))
        g = np.random.default_rng(SEED + 7)
        s1, s2 = [], []
        for _ in range(B_SIM):
            noisy = mapped + g.normal(0, sigma, len(mapped))
            s1.append(sp(paH1[m], noisy))
            s2.append(sp(paH1[m] + g.normal(0, sigmaH1, m.sum()), noisy))
        s1 = np.array([v for v in s1 if np.isfinite(v)])
        s2 = np.array([v for v in s2 if np.isfinite(v)])
        obs_m = sp(paH1[m], tgt[m])
        b1 = (float(np.percentile(s1, 2.5)), float(np.percentile(s1, 97.5)))
        b2 = (float(np.percentile(s2, 2.5)), float(np.percentile(s2, 97.5)))
        rows.append(dict(
            sample=s, obs_rho=obs_m, band_lo=b1[0], band_hi=b1[1],
            sim_med=float(np.median(s1)),
            in_band=bool(b1[0] <= obs_m <= b1[1]),
            direction=("inside" if b1[0] <= obs_m <= b1[1]
                       else ("below" if obs_m < b1[0] else "above")),
            sigma=sigma, r_full=r_full, r_full_old=float(old.loc[s, "r_full"]),
            band_lo_2side=b2[0], band_hi_2side=b2[1],
            in_band_2side=bool(b2[0] <= obs_m <= b2[1]),
            sigma_H1=sigmaH1,
            pa5_median=float(s4.loc[s, "med_pa5"]) if s in s4.index else np.nan,
            is_pluri=bool(s in PLURI)))
    S = pd.DataFrame(rows)
    S.to_csv(OUT, sep="\t", index=False)
    print(S.round(4).to_string(index=False))
    print("\n样本数", len(S), "→", OUT)
    for s in ("H9", "WTC-11", "H1"):
        if s in S["sample"].values or s in rel.index:
            r = S[S["sample"] == s]
            if len(r):
                r = r.iloc[0]
                print(f"{s}: obs {r.obs_rho:.3f} band [{r.band_lo:.3f}, {r.band_hi:.3f}] "
                      f"excess {r.obs_rho - r.band_hi:+.4f}")


if __name__ == "__main__":
    main()
