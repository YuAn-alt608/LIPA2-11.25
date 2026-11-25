#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
311 — 阈值敏感性（ROUND13 P0-3）+ 样本两两一致性（P0-4）。

P0-3
----
阈值 0.5 / 0.6 / 0.667 / 0.7 / 0.75 / 0.8 / 0.9，对 15 个有单倍型轨道的样本：
  点估计可读数、区间完全高于阈值、borderline（区间跨阈值）、不可读。
H1 是参考样本，原表 r_full 的区间是占位（1.0/1.0），这里用 23 染色体块 bootstrap
重算 r_s_nz -> r_full 的区间（B=1000，每个重复内重新去双零）。

P0-4
----
11 样本面板所有 55 对：ρ_obs（L1PA2 1,032 拷贝的 pa5），有 r_full 的样本对再给
衰减校正 ρ_corr = ρ_obs / sqrt(r_full,i * r_full,j)。区间用 23 染色体块 bootstrap，
每个重复内重算 r_s_nz -> r_full -> ρ_obs -> ρ_corr（B=1000）。
GM12878 没有单倍型轨道，与它相关的对只报 ρ_obs。

输出：t13/t13_threshold_sensitivity.tsv
      t13/t13_pairwise_rho.tsv
      t13/t13_pairwise_summary.tsv
      t13/t13_h1_rfull_ci.tsv
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from itertools import combinations

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
OUT = os.path.join(RES, "t13")
T9B = os.path.join(RES, "t9b")
T10 = RES
B = 1000
SEED = 20261124
PLURI = ["H1", "H9", "WTC-11", "GM27730"]
NONPL = ["K562", "HepG2", "Hap1", "Hek293T", "MCF-7", "HG002", "GM12878"]
S11 = PLURI + NONPL
READ_NONPL = ["Hek293T", "MCF-7", "HepG2", "K562"]
THR = [0.5, 0.6, 0.667, 0.7, 0.75, 0.8, 0.9]


def sb(r):
    return 2.0 * r / (1.0 + r) if np.isfinite(r) and r > -1 else np.nan


def sp(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def blocks(chrom):
    c = pd.Categorical(np.asarray(chrom, str))
    return [np.where(c.codes == j)[0] for j in range(len(c.categories))]


def r_nz(a1, a2):
    a1 = np.asarray(a1, float); a2 = np.asarray(a2, float)
    ok = np.isfinite(a1) & np.isfinite(a2)
    nz = ok & ~((a1 == 0) & (a2 == 0))
    return sp(a1[nz], a2[nz])


def main():
    os.makedirs(OUT, exist_ok=True)
    rel = pd.read_csv(os.path.join(T10, "t10_reliability_fig.tsv"), sep="\t")
    hap = pd.read_csv(os.path.join(T9B, "t9b_pa5_hap_L1PA2_fixed.tsv.gz"), sep="\t")
    pa = pd.read_csv(os.path.join(RES, "t1d_pa5_by_copy_11.tsv.gz"), sep="\t")
    pa = pa[pa.grp == "L1PA2"].reset_index(drop=True)
    assert len(pa) == 1032 and len(hap) == 1032, (len(pa), len(hap))
    # 行序断言：hap 与 pa 必须逐行同序（都源自 t9b 修复表）
    m = pa.merge(hap, on="unit_id", how="inner", validate="one_to_one")
    assert len(m) == 1032, "hap 与 pa 的 unit_id 不一一对应"
    hap = m
    chrom = hap.unit_id.str.split("_").str[0].to_numpy()

    # ---------- H1 的 r_full 区间 ----------
    a1 = hap["H1.hap1"].to_numpy(float); a2 = hap["H1.hap2"].to_numpy(float)
    r0 = r_nz(a1, a2)
    bb = blocks(chrom)
    rng = np.random.default_rng(SEED)
    vals = [sb(r_nz(a1[ix], a2[ix])) for ix in
            (np.concatenate([bb[j] for j in rng.integers(0, len(bb), len(bb))]) for _ in range(B))]
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    H1 = pd.DataFrame([dict(sample="H1", r_s_nz=r0, r_full=sb(r0), r_full_lo=lo, r_full_hi=hi,
                            n=1032, B=B, note="23-block bootstrap, 每个重复内重新去双零")])
    H1.to_csv(os.path.join(OUT, "t13_h1_rfull_ci.tsv"), sep="\t", index=False)
    print(f"H1 r_s_nz={r0:.4f} r_full={sb(r0):.4f} [{lo:.4f}, {hi:.4f}]", flush=True)

    # ---------- P0-3 阈值表 ----------
    rel = rel[rel["sample"] != "Hap1"].copy()          # 无效，不参与计数
    rel.loc[rel["sample"] == "H1", "r_full_lo"] = lo
    rel.loc[rel["sample"] == "H1", "r_full_hi"] = hi
    rows = []
    for t in THR:
        v = rel.dropna(subset=["r_full"])
        point = int((v.r_full >= t).sum())
        above = int(((v.r_full_lo >= t) & (v.r_full_hi >= t)).sum())
        border = int(((v.r_full_lo < t) & (v.r_full_hi >= t)).sum())
        below = int((v.r_full_hi < t).sum())
        names_above = ", ".join(sorted(v.loc[v.r_full_lo >= t, "sample"]))
        names_border = ", ".join(sorted(v.loc[(v.r_full_lo < t) & (v.r_full_hi >= t), "sample"]))
        rows.append(dict(threshold=t, n_assessed=len(v), point_readable=point,
                         ci_fully_above=above, borderline=border, unreadable=below,
                         ci_above_samples=names_above, borderline_samples=names_border))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "t13_threshold_sensitivity.tsv"), sep="\t", index=False)

    # ---------- P0-4 两两一致性 ----------
    M = pa[S11].to_numpy(float)
    relmap = rel.set_index("sample")
    hapcols = {s: (hap[f"{s}.hap1"].to_numpy(float), hap[f"{s}.hap2"].to_numpy(float))
               for s in S11 if f"{s}.hap1" in hap.columns}
    rfull0 = {s: (float(relmap.loc[s, "r_full"]) if s in relmap.index else np.nan) for s in S11}

    rng = np.random.default_rng(SEED + 1)
    picks = [np.concatenate([bb[j] for j in rng.integers(0, len(bb), len(bb))]) for _ in range(B)]

    # 每个重复内先重算各样本 r_full（同一次重采样索引，保持样本间配对）
    rf_boot = {s: np.full(B, np.nan) for s in hapcols}
    for s, (x, y) in hapcols.items():
        for b, ix in enumerate(picks):
            rf_boot[s][b] = sb(r_nz(x[ix], y[ix]))

    rows = []
    for s1, s2 in combinations(S11, 2):
        x, y = M[:, S11.index(s1)], M[:, S11.index(s2)]
        rho0 = sp(x, y)
        rf1, rf2 = rfull0[s1], rfull0[s2]
        corr0 = rho0 / np.sqrt(rf1 * rf2) if np.isfinite(rf1) and np.isfinite(rf2) and rf1 > 0 and rf2 > 0 else np.nan
        # ρ_obs 的块 bootstrap
        rb = np.array([sp(x[ix], y[ix]) for ix in picks])
        lo_o, hi_o = np.nanpercentile(rb, [2.5, 97.5])
        # ρ_corr 的块 bootstrap：同一重复内用该重复的 r_full
        if np.isfinite(corr0):
            cb = rb / np.sqrt(rf_boot[s1] * rf_boot[s2])
            lo_c, hi_c = np.nanpercentile(cb, [2.5, 97.5])
        else:
            lo_c = hi_c = np.nan
        g1 = "pluripotent" if s1 in PLURI else "other"
        g2 = "pluripotent" if s2 in PLURI else "other"
        group = ("pluripotent x pluripotent" if g1 == g2 == "pluripotent" else
                 "other x other" if g1 == g2 == "other" else "pluripotent x other")
        rows.append(dict(sample_a=s1, sample_b=s2, group=group,
                         rho_obs=rho0, rho_obs_lo=lo_o, rho_obs_hi=hi_o,
                         r_full_a=rf1, r_full_b=rf2,
                         rho_corr=corr0, rho_corr_lo=lo_c, rho_corr_hi=hi_c,
                         n_copy=1032))
    PR = pd.DataFrame(rows)
    PR.to_csv(os.path.join(OUT, "t13_pairwise_rho.tsv"), sep="\t", index=False)

    # 分组汇总（按稿件 P0-4 的三组；多能×可读非多能只取 4 个样本）
    def summ(sub, name):
        v = sub.rho_corr.dropna()
        return dict(group=name, n_pairs=len(sub), n_with_corr=len(v),
                    min=v.min(), median=v.median(), max=v.max(),
                    obs_min=sub.rho_obs.min(), obs_median=sub.rho_obs.median(),
                    obs_max=sub.rho_obs.max())
    sub_p = PR[(PR.sample_a.isin(PLURI)) & (PR.sample_b.isin(PLURI))]
    sub_pn = PR[((PR.sample_a.isin(PLURI)) & (PR.sample_b.isin(READ_NONPL))) |
                ((PR.sample_b.isin(PLURI)) & (PR.sample_a.isin(READ_NONPL)))]
    sub_n = PR[(PR.sample_a.isin(READ_NONPL)) & (PR.sample_b.isin(READ_NONPL))]
    SS = pd.DataFrame([summ(sub_p, "pluripotent internal"),
                       summ(sub_pn, "pluripotent x readable non-pluripotent"),
                       summ(sub_n, "readable non-pluripotent internal")])
    SS.to_csv(os.path.join(OUT, "t13_pairwise_summary.tsv"), sep="\t", index=False)

    # 多能内部不含 H1 的 3 对
    sub3 = sub_p[(sub_p.sample_a != "H1") & (sub_p.sample_b != "H1")]
    print(SS.to_string(), flush=True)
    print("pluri pairs excluding H1: corr", round(sub3.rho_corr.min(), 3),
          round(sub3.rho_corr.median(), 3), round(sub3.rho_corr.max(), 3), flush=True)

    L = ["# 阈值敏感性与两两一致性（311 / P0-3、P0-4）", "",
         f"H1 的 r_full 区间由 23 染色体块 bootstrap 得到（B={B}，每个重复内重新去双零）："
         f"r_s_nz = {r0:.3f}，r_full = {sb(r0):.3f} [{lo:.3f}, {hi:.3f}]。", "",
         "## 阈值敏感性（15 个有单倍型轨道的样本，Hap1 无效已排除）", "",
         "| 阈值 | 点估计可读 | 区间完全高于 | borderline | 不可读 | 区间高于的样本 |", "|---|---|---|---|---|---|"]
    for _, r in T.iterrows():
        L.append(f"| {r.threshold} | {r.point_readable} | {r.ci_fully_above} | {r.borderline} | "
                 f"{r.unreadable} | {r.ci_above_samples} |")
    L += ["", "## 两两一致性（按组）", "",
          "| 组 | 对数 | 有校正值 | ρ_corr 最小 | 中位 | 最大 | ρ_obs 最小–最大 |", "|---|---|---|---|---|---|---|"]
    for _, r in SS.iterrows():
        L.append(f"| {r.group} | {r.n_pairs} | {r.n_with_corr} | {r['min']:.3f} | {r['median']:.3f} | "
                 f"{r['max']:.3f} | {r.obs_min:.2f}–{r.obs_max:.2f} |")
    open(os.path.join(OUT, "t13_threshold_pairwise.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
