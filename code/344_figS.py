#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
344 — 补充图 S1–S5。

S1 c5 全部 11 个阈值（0.70/0.75 在预设范围外，空心）：上=名次，下=perm_p（对数）
S2 13 个样本的 ρ_obs 与噪声带
S3 定相诊断：双零比例 / 池化-单倍型相关 与 均值-池化相关
S4 等位基因：946 拷贝 ΔM vs ΔA；|ΔA| 五分位
S5 31 个「样本×家族」格子的 r_full vs pa5
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
fc = importlib.import_module("340_fig_common")


def ax_mm(fig, x, y, w, h, W, H):
    return fig.add_axes([x / W, y / H, w / W, h / H])


# --------------------------------------------------------------- S1
def fig_s1():
    print("Fig. S1")
    c5 = fc.load("t1_did_perm4_c5.tsv").sort_values("threshold")
    W, H = 174.0, 96.0
    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))
    ax = ax_mm(fig, 12.0, 52.0, 162.0, 36.0, W, H)
    x = c5.threshold.values
    xin = c5.in_prespecified_range.astype(bool).values
    for xx in x[~xin]:
        ax.axvspan(xx - 0.014, xx + 0.014, color="#F0F0F0", zorder=0)
    for s_, mk in (("H1", "o"), ("H9", "s"), ("WTC11", "^"), ("GM27730", "D")):
        lab = "WTC-11" if s_ == "WTC11" else s_
        v = c5[f"rank_{s_}"].values
        ax.plot(x[xin], v[xin], marker=mk, ms=2.8, color=fc.BLUE, lw=0.6, mfc=fc.BLUE,
                mew=0.6, label=lab)
        ax.plot(x[~xin], v[~xin], marker=mk, ms=2.8, color=fc.BLUE, lw=0.6, mfc="white",
                mew=0.7, ls=":")
    ax.set_yscale("log")
    ax.set_ylim(46, 0.85)
    ax.set_yticks([1, 2, 4, 8, 16, 32])
    ax.set_yticklabels(["1", "2", "4", "8", "16", "32"])
    ax.set_xlim(0.23, 0.78)
    ax.set_xticks([])
    ax.set_ylabel("DiD rank (log)")
    fc.clean(ax)
    ax.legend(loc="lower left", frameon=False, ncol=4, handletextpad=0.3, columnspacing=1.6,
              fontsize=8, borderaxespad=0.2)
    fig.text(0.50, 0.505, "grey band = outside the pre-specified range", fontsize=8,
             ha="center", va="top", color=fc.GREY)

    ax2 = ax_mm(fig, 12.0, 16.0, 162.0, 28.0, W, H)
    ax2.plot(x[xin], c5.perm_p.values[xin], marker="o", ms=2.8, color=fc.GREY, lw=0.7)
    ax2.plot(x[~xin], c5.perm_p.values[~xin], marker="o", ms=2.8, mfc="white", mec=fc.GREY,
             lw=0.7, ls=":")
    ax2.set_yscale("log")
    ax2.set_xlim(0.23, 0.78)
    ax2.set_ylim(5e-6, 2e-4)
    ax2.set_yticks([1e-5, 1e-4])
    ax2.set_xticks(list(x))
    ax2.set_xticklabels([f"{v:.2f}" for v in x], rotation=45, ha="right", fontsize=8)
    ax2.set_xlabel("c5 cut-off", labelpad=1)
    ax2.set_ylabel("permutation p")
    fc.clean(ax2)
    ax2.text(0.775, 6.0e-6, f"C(41,4) = {int(c5.n_comb.iloc[0]):,}", fontsize=8,
             color=fc.GREY, ha="right")

    c5.to_csv(os.path.join(fc.FIGDIR, "figS1_data_c5_all.tsv"), sep="\t", index=False)
    fc.dump_numbers("figs1", {
        "n_thresholds": int(len(c5)),
        "prespecified": [float(v) for v in x[xin]],
        "outside": [float(v) for v in x[~xin]],
        "perm_p_prespecified": [float(c5.perm_p.values[xin].min()), float(c5.perm_p.values[xin].max())],
        "n_comb": int(c5.n_comb.iloc[0]),
        "top4_pluripotent_prespecified": int(c5.top4_is_group.astype(bool).values[xin].sum()),
    })
    fc.save(fig, "FigS1", 174)


# --------------------------------------------------------------- S2
def fig_s2():
    print("Fig. S2")
    b = fc.load("t10_sim_band_cov.tsv").sort_values("obs_rho", ascending=False).reset_index(drop=True)
    # 40 个样本逐行标注：行距需 >= 4.4 mm，否则 8 pt 的 y 标签互相压字
    W, H = 174.0, 200.0
    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))
    ax = ax_mm(fig, 15.0, 12.0, 157.0, 176.0, W, H)
    yy = np.arange(len(b))[::-1]
    for y, (_, r) in zip(yy, b.iterrows()):
        c = fc.colour(r["sample"])
        ax.add_patch(Rectangle((r.band_lo, y - 0.26), r.band_hi - r.band_lo, 0.52,
                               facecolor="#E6E6E6", edgecolor="none", zorder=0))
        ax.plot([r.obs_rho], [y], marker="o", ms=3.0, color=c, zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(fc.disp(list(b["sample"])), fontsize=8)
    ax.tick_params(axis="y", length=0, pad=1)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.7, len(b) - 0.2)
    ax.set_xlabel("ρ_obs against H1", labelpad=1)
    fc.clean(ax)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=fc.GREY, ms=3.0, label="ρ_obs"),
                       Rectangle((0, 0), 1, 1, facecolor="#E6E6E6", edgecolor="none",
                                 label="noise band")],
              loc="lower right", frameon=False, handletextpad=0.3, labelspacing=0.2,
              fontsize=8, borderaxespad=0.1)
    h9 = b[b["sample"] == "H9"].iloc[0]
    ax.annotate(f"+{h9.obs_rho - h9.band_hi:.3f}", (h9.obs_rho, len(b) - 1 + 0.34),
                fontsize=8, ha="center", va="bottom", color=fc.GREY)
    ax.text(0.02, 0.97, f"H9 self-comparison {h9.obs_rho - h9.band_hi:.3f} above its band;\nreadable non-pluripotent samples below theirs",
            transform=ax.transAxes, fontsize=8, va="top")
    b.to_csv(os.path.join(fc.FIGDIR, "figS2_data_band.tsv"), sep="\t", index=False)
    fc.dump_numbers("figs2", {
        "h9_excess": float(h9.obs_rho - h9.band_hi),
        "n_below": int((b.direction == "below").sum()),
        "n_inside": int((b.direction == "inside").sum()),
        "n_samples": int(len(b)),
    })
    fc.save(fig, "FigS4", 174)


# --------------------------------------------------------------- S3
def fig_s3():
    print("Fig. S3")
    rel = fc.rel41()
    sc = fc.load("t15/t15_selfconsistency.tsv")
    # (a) 画的是有效性规则第一条检验的量：两条单倍型都没覆盖的拷贝占比
    _cv = fc.load("t15/t15_reliability_cov_41.tsv")[["sample", "frac_nocov"]]
    rel = rel.merge(_cv, on="sample", how="left")
    rel["frac_both0"] = pd.to_numeric(rel["frac_both0"], errors="coerce")
    rel["frac_nocov"] = pd.to_numeric(rel["frac_nocov"], errors="coerce")
    rel["rho_m_h1"] = pd.to_numeric(rel["rho_m_h1"], errors="coerce")
    rel["rho_m_h2"] = pd.to_numeric(rel["rho_m_h2"], errors="coerce")
    order = fc.sample_order_by_rfull(rel)
    rel = rel.set_index("sample").loc[order].reset_index()
    rel["y"] = np.arange(len(rel))[::-1]

    W, H = 174.0, 168.0
    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))
    axa = ax_mm(fig, 15.0, 18.0, 50.0, 136.0, W, H)
    for _, r in rel.iterrows():
        c = fc.colour(r["sample"])
        axa.plot([0, r.frac_nocov], [r.y, r.y], color=fc.LGREY, lw=0.5, zorder=1)
        axa.plot([r.frac_nocov], [r.y], marker="o", ms=2.4, color=c, zorder=3)
    axa.axvline(0.5, color=fc.GREY, lw=0.6, ls=(0, (3, 2)))
    axa.text(0.507, len(rel) - 0.55, "0.5", fontsize=8, color=fc.GREY, va="top")
    axa.set_yticks(rel.y)
    axa.set_yticklabels(fc.disp(list(rel["sample"])), fontsize=8)
    # 有 11 个样本的 frac_nocov 只有 0.005–0.01，圆点几乎贴在轴线上，
    # pad 太小会让圆点和左侧标签看起来连在一起
    axa.tick_params(axis="y", length=0, pad=4)
    axa.set_xlim(0, 1.0)
    axa.set_ylim(-0.8, len(rel) - 0.2)
    axa.set_xticks([0, 0.5, 1.0])
    axa.set_xlabel("fraction with neither haplotype covered", labelpad=1)
    fc.clean(axa)
    fc.panel(axa, "a", dx=0.0, dy=1.005)

    axb = ax_mm(fig, 82.0, 18.0, 92.0, 136.0, W, H)
    sc = sc.set_index("sample")
    for _, r in rel.iterrows():
        c = fc.colour(r["sample"])
        phased = bool(r.phased)
        if np.isfinite(r.rho_m_h1) and np.isfinite(r.rho_m_h2):
            # 未通过定相有效性判据的两个样本画空心（图例不再用 ×，因为图上没有 ×）
            axb.plot([min(r.rho_m_h1, r.rho_m_h2)], [r.y], marker="o", ms=2.3,
                     mfc=c if phased else "white", mec=c, mew=0.6 if phased else 0.9, zorder=3)
        if phased and r["sample"] in sc.index and np.isfinite(sc.loc[r["sample"], "rho_mean_vs_pooled"]):
            axb.plot([sc.loc[r["sample"], "rho_mean_vs_pooled"]], [r.y], marker="s", ms=2.0,
                     mfc="white", mec=c, mew=0.6, zorder=3)
    axb.axvline(0.2, color=fc.GREY, lw=0.6, ls=(0, (3, 2)))
    axb.text(0.212, len(rel) * 0.46, "0.2", fontsize=8, color=fc.GREY, va="center", ha="left")
    axb.set_yticks([])
    axb.set_xlim(-0.12, 1.02)
    axb.set_ylim(-0.8, len(rel) - 0.2)
    axb.set_xticks([0, 0.5, 1.0])
    axb.set_xlabel("correlation with pooled track", labelpad=1)
    fc.clean(axb, left=False)
    axb.legend(handles=[Line2D([], [], marker="o", ls="", color=fc.GREY, ms=2.6,
                               label="min(ρ pooled–hap)"),
                        Line2D([], [], marker="s", ls="", mfc="white", mec=fc.GREY, mew=0.6,
                               ms=2.4, label="mean of haps vs pooled"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec=fc.GREY, mew=0.9,
                               ms=2.6, label="not assessable")],
               loc="upper left", frameon=False, fontsize=8, handletextpad=0.35,
               labelspacing=0.2, borderaxespad=0.1)
    fc.panel(axb, "b", dx=-0.035, dy=1.005)


    rel.to_csv(os.path.join(fc.FIGDIR, "figS3_data_phasing.tsv"), sep="\t", index=False)
    v = rel[np.isfinite(rel.rho_m_h1) & np.isfinite(rel.rho_m_h2)]
    ph = set(rel.loc[rel.phased.astype(bool), "sample"])
    fc.dump_numbers("figs3", {
        "n_samples": int(len(rel)), "n_phased": int(rel.phased.sum()),
        "frac_both0_max": float(rel.frac_both0.max()),
        "hap1_frac_nocov": float(rel.set_index("sample").loc["Hap1", "frac_nocov"]),
        "hap1_frac_both0": float(rel.set_index("sample").loc["Hap1", "frac_both0"]),
        "gm12878_rho_h2": float(rel.set_index("sample").loc["GM12878", "rho_m_h2"]),
        "gm12878_nz_hap2": float(rel.set_index("sample").loc["GM12878", "nz_hap2"]),
        "selfcons_range": [float(sc.loc[sc.index.isin(ph), "rho_mean_vs_pooled"].min()),
                           float(sc.loc[sc.index.isin(ph), "rho_mean_vs_pooled"].max())],
        "selfcons_median": float(sc.loc[sc.index.isin(ph), "rho_mean_vs_pooled"].median()),
    })
    fc.save(fig, "FigS2", 174)


# --------------------------------------------------------------- S4
def fig_s4():
    print("Fig. S4")
    d = fc.load("t9i/t9i_dA_by_copy.tsv.gz")
    m = d.dA_hap_H1.notna() & d.dM_hap_H1.notna() & (d.nCpg_hap1 > 0) & (d.nCpg_hap2 > 0)
    s = d[m].copy()
    W, H = 174.0, 88.0
    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))
    axa = ax_mm(fig, 9.0, 16.0, 80.0, 62.0, W, H)
    axa.plot(s.dA_hap_H1, s.dM_hap_H1, ls="", marker="o", ms=1.2, mfc=fc.GREY, mec="none", alpha=0.45)
    axa.axhline(0, color=fc.LGREY, lw=0.5)
    axa.axvline(0, color=fc.LGREY, lw=0.5)
    axa.set_xlim(-32, 33)
    axa.set_ylim(-0.55, 0.42)
    axa.set_xticks([-30, -15, 0, 15, 30])
    axa.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])
    axa.set_xlabel("ΔA (pa5, hap1 − hap2, pp)", labelpad=1)
    axa.set_ylabel("ΔM (5mC, hap1 − hap2)", labelpad=1)
    fc.clean(axa)
    from scipy.stats import spearmanr
    rho = float(spearmanr(s.dA_hap_H1, s.dM_hap_H1).statistic)
    axa.text(0.03, 0.97, f"ρ = {fc.num(rho)} (95% CI \u22120.33 to \u22120.22)\nn = {len(s)}",
             transform=axa.transAxes, fontsize=8, va="top")
    fc.panel(axa, "a", dx=0.0, dy=1.01)

    axb = ax_mm(fig, 101.0, 16.0, 73.0, 62.0, W, H)
    q = pd.qcut(s.dA_hap_H1.abs(), 5, labels=False)
    vals = [float(spearmanr(g.dA_hap_H1, g.dM_hap_H1).statistic) for _, g in s.groupby(q)]
    axb.axhline(0, color=fc.LGREY, lw=0.5)
    axb.bar(np.arange(5), vals, width=0.62, color=fc.BLUE, edgecolor="none", zorder=2)
    axb.set_xticks(np.arange(5))
    axb.set_xticklabels(["Q1\nlow", "Q2", "Q3", "Q4", "Q5\nhigh"], fontsize=8)
    axb.set_xlim(-0.6, 4.6)
    axb.set_ylim(-0.62, 0.06)
    axb.set_yticks([-0.6, -0.4, -0.2, 0])
    axb.set_xlabel("quintile of |ΔA|", labelpad=1)
    axb.set_ylabel("ρ(ΔA, ΔM)", labelpad=1)
    fc.clean(axb)
    fc.panel(axb, "b", dx=-0.05, dy=1.01)
    s.to_csv(os.path.join(fc.FIGDIR, "figS4_data_allele.tsv"), sep="\t", index=False)
    fc.dump_numbers("figs4", {"n": int(len(s)), "rho": rho, "quintiles": vals,
                              "opposite_sign_frac": float((np.sign(s.dA_hap_H1) != np.sign(s.dM_hap_H1)).mean()),
                              "sd_dA": float(s.dA_hap_H1.std()), "sd_dM": float(s.dM_hap_H1.std())})
    fc.save(fig, "FigS5", 174)


# --------------------------------------------------------------- S5
def fig_s5():
    print("Fig. S5")
    d = fc.load("t10/t10_readability_map.tsv")
    d = d[(d.phase_valid == True) & d.r_s_nz.notna()].copy()
    print("   cells:", len(d), d.family.value_counts().to_dict())
    assert len(d) == 31, len(d)
    W, H = 174.0, 118.0
    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))
    ax = ax_mm(fig, 9.0, 16.0, 165.0, 92.0, W, H)
    famcol = {"L1PA2": fc.BLUE, "L1HS": fc.GREEN, "LTR7": fc.ORANGE, "AluS": fc.PURPLE}
    for f, g in d.groupby("family"):
        ax.plot(g.pa5_median, g.r_full, ls="", marker="o", ms=3.2,
                color=famcol.get(f, fc.GREY), label=f)
    ax.set_xscale("log")
    ax.set_xlim(0.8, 25)
    ax.set_ylim(-0.03, 1.0)
    ax.set_xticks([1, 2, 5, 10, 20])
    ax.set_xticklabels(["1", "2", "5", "10", "20"])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("pa5 median (%), log scale", labelpad=1)
    ax.set_ylabel("r_full", labelpad=1)
    fc.clean(ax)
    from scipy.stats import spearmanr
    rho31 = float(spearmanr(d.pa5_median, d.r_full).statistic)
    sub = d[d.family.isin(["L1HS", "LTR7", "AluS"])]
    rho22 = float(spearmanr(sub.pa5_median, sub.r_full).statistic)
    ax.text(0.02, 0.05, f"31 cells: +{rho31:.2f}\nL1HS/LTR7/AluS ({len(sub)}): +{rho22:.2f}",
            transform=ax.transAxes, fontsize=8, va="bottom")
    ax.legend(loc="upper left", frameon=False, ncol=4, handletextpad=0.2, columnspacing=1.4,
              labelspacing=0.25, fontsize=8)
    d.to_csv(os.path.join(fc.FIGDIR, "figS5_data_cells.tsv"), sep="\t", index=False)
    fc.dump_numbers("figs5", {"n_cells": int(len(d)), "by_family": d.family.value_counts().to_dict(),
                              "rho31": rho31, "n_sub": int(len(sub)), "rho_sub": rho22})
    fc.save(fig, "FigS3", 174)


if __name__ == "__main__":
    fig_s1(); fig_s2(); fig_s3(); fig_s4(); fig_s5()
