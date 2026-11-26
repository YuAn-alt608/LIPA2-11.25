#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
341 — Fig. 1（174 × 132 mm；六个面板，按阅读顺序 a–f）

字母顺序 = 正文首次引用顺序：
(a) c5 九个预设阈值的 DiD 名次曲线（Results 段 1）
(b) pa5 口径的名次（九样本面板，四点标样本名）（Results 段 1）
(c) 逐样本 pa5 中位数（5′ 1 kb）（Results 段 2）
(d) 四家族 raw 与 DiD（5′ 1 kb）（Results 段 2）
(e) 逐样本 L1PA2−L1PA4–7 差值 + 留一带（Results 段 3）
(f) 400 bp 窗的对照家族森林图（Results 段 4）

r22 改动：(c)(d) 与 (f) 换位——原 c=400 bp 森林、d=median、f=四家族，正文却按
d→f→e→c 引用；(b) 补四个样本名；森林图内的 positive/near-zero control 文字移到图注
（面板变窄后放不下）。

版面按毫米用 add_axes 定位，避免 GridSpec 在窄栏里挤压相邻面板。
所有数值来自 03_results；跑前做 3 拷贝 bigWig 抽查。
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
fc = importlib.import_module("340_fig_common")

RES = fc.RES
W, H = 174.0, 132.0


def ax_mm(fig, x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


def main():
    print("Fig. 1")
    fc.spot_check("fig1")

    c5 = fc.load("t1_did_perm4_c5.tsv")
    c5 = c5[c5.in_prespecified_range.astype(bool)].sort_values("threshold")
    pa5p = fc.load("t1_did_perm4_pa5.tsv").iloc[0]
    gc = fc.load("t13/t13_group_contrast_sample_first.tsv")
    ps = fc.load("t13/t13_per_sample_did.tsv")

    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))

    # ---------------- (a) 左上：名次曲线 ----------------
    ax = ax_mm(fig, 13.5, 78.0, 49.5, 46.0)
    x = c5.threshold.values
    for s, mk, filled in (("H1", "o", True), ("H9", "s", True),
                          ("WTC11", "^", True), ("GM27730", "D", False)):
        lab = "WTC-11" if s == "WTC11" else s
        ax.plot(x, c5[f"rank_{s}"].values, marker=mk, ms=2.6, color=fc.BLUE, lw=0.6,
                mfc=fc.BLUE if filled else "white", mew=0.7, label=lab, clip_on=False)
    ax.set_xlim(0.235, 0.685)
    ax.set_yticks([1, 2, 3, 4, 5, 7, 9])
    ax.set_xticks([0.25, 0.35, 0.45, 0.55, 0.65])
    ax.set_xlabel("c5 cut-off", labelpad=1)
    ax.set_ylabel("DiD rank (1 = highest)", labelpad=2)
    fc.clean(ax)
    for _, r in c5[~c5.top4_is_group.astype(bool)].iterrows():
        ax.text(r.threshold, 0.30, str(r.intruder), rotation=90, fontsize=8,
                ha="center", va="bottom", color="#000000")
    # 名次 9 的两个点（0.55、0.65）与 GM27730 的下降臂会切过图例，
    # 故把纵轴下沿放宽到 13.6 腾出底部空白：注释在上、图例在下
    ax.set_ylim(13.6, 0.15)
    ax.text(0.02, 0.175, f"top-4 = pluripotent at {int(c5.top4_is_group.sum())}/{len(c5)} cut-offs",
            transform=ax.transAxes, fontsize=8, va="bottom", ha="left")
    ax.legend(loc="lower left", bbox_to_anchor=(0.005, 0.005), frameon=False, ncol=2,
              handletextpad=0.2, columnspacing=0.8, labelspacing=0.15, borderpad=0.1,
              fontsize=8)
    fc.panel(ax, "a", dx=-0.26, dy=1.01)

    # ---------------- (b) 中上：pa5 名次（四点标名） ----------------
    axp = ax_mm(fig, 68.0, 78.0, 26.0, 46.0)
    nine = fc.load("t1d_pa5_by_copy_9.tsv.gz")
    nine = nine[nine.grp == "L1PA2"]
    cols9 = ["H1", "WTC-11", "H9", "GM27730", "Hek293T", "MCF-7", "HepG2", "Hap1", "K562"]
    med9 = nine[cols9].median().sort_values(ascending=False)
    PLURI = {"H1", "H9", "WTC-11", "GM27730"}
    for rk, (name, val) in enumerate(med9.items(), start=1):
        is_p = name in PLURI
        axp.plot([0.22], [rk], marker="o", ms=3.2,
                 color=fc.BLUE if is_p else fc.GREY, clip_on=False)
        if is_p:
            axp.text(0.40, rk, fc.disp(name), fontsize=8, va="center", ha="left", color="#000000")
    axp.set_xlim(0, 2.6)
    axp.set_ylim(13.6, 0.15)
    axp.set_yticks([1, 2, 3, 4])
    axp.set_xticks([])
    axp.tick_params(axis="y", pad=1.5, labelsize=8)
    fc.clean(axp)
    axp.set_xlabel("pa5\n(n = 9)", labelpad=1)
    # 精确 p 不画在图上（无空白可用，放图注）：1/C(9,4) = 1/126 = 0.0079
    fc.panel(axp, "b", dx=-0.07, dy=1.02)

    # ---------------- (c) 右上：逐样本中位数 ----------------
    order = list(ps[ps.group == "pluripotent"].sort_values("diff", ascending=False)["sample"]) + \
            list(ps[ps.group != "pluripotent"].sort_values("diff", ascending=False)["sample"])
    ny = len(order)
    ax1 = ax_mm(fig, 118.0, 78.0, 42.0, 46.0)
    for i, smp in enumerate(order):
        r = ps[ps["sample"] == smp].iloc[0]
        ax1.plot([0, r.l1pa2_median], [i, i], color=fc.LGREY, lw=0.45, zorder=1)
        ax1.plot([r.l1pa2_median], [i], marker="o", ms=3.0, color=fc.colour(smp), zorder=3)
    gpl = ps[ps.group == "pluripotent"].l1pa2_median.mean()
    got = ps[ps.group != "pluripotent"].l1pa2_median.mean()
    ax1.plot([gpl, gpl], [-0.55, 3.5], color=fc.BLUE, lw=0.7, ls=(0, (2.5, 1.5)), zorder=2)
    ax1.plot([got, got], [3.5, ny - 0.45], color=fc.GREY, lw=0.7, ls=(0, (2.5, 1.5)), zorder=2)
    ax1.text(gpl + 0.9, 3.35, f"{gpl:.1f}", fontsize=8, color=fc.BLUE, ha="left", va="bottom")
    ax1.text(8.6, 10.3, f"{got:.1f}", fontsize=8, color=fc.GREY, ha="left", va="center")
    ax1.axhline(3.5, color=fc.LGREY, lw=0.5, zorder=0)
    ax1.set_xlim(0, 22)
    ax1.set_ylim(ny - 0.45, -0.55)
    ax1.set_xticks([0, 5, 10, 15, 20])
    ax1.set_yticks(np.arange(ny))
    ax1.set_yticklabels([fc.disp(x) for x in order], fontsize=8)
    ax1.tick_params(axis="y", length=0, pad=2)
    ax1.set_xlabel("pa5 median (%)", labelpad=2)
    fc.clean(ax1)
    ax1.set_title("median", fontsize=8, pad=3)
    fc.panel(ax1, "c", dx=-0.40, dy=1.10)

    # ---------------- (d) 左下：四家族森林 ----------------
    ax2 = ax_mm(fig, 29.0, 16.0, 40.0, 46.0)
    fams = ["L1HS", "L1PA2", "L1PA3", "L1PA4-7"]
    yp = np.arange(len(fams))[::-1]
    for y, f in zip(yp, fams):
        r = gc[(gc.window == "5p1k") & (gc.family == f)].iloc[0]
        ax2.errorbar(r.raw_diff, y + 0.17, xerr=[[r.raw_diff - r.raw_lo], [r.raw_hi - r.raw_diff]],
                     fmt="o", ms=3.0, mfc="white", mec=fc.GREY, ecolor=fc.GREY, elinewidth=0.6,
                     capsize=1.2, zorder=3)
        if np.isfinite(r.did):
            ax2.errorbar(r.did, y - 0.17, xerr=[[r.did - r.did_lo], [r.did_hi - r.did]],
                         fmt="o", ms=3.2, color=fc.BLUE, ecolor=fc.BLUE, elinewidth=0.6,
                         capsize=1.2, zorder=3)
    ax2.axvline(0, color=fc.LGREY, lw=0.5, zorder=1)
    ax2.set_yticks(yp)
    ax2.set_yticklabels(fams, fontsize=8)
    ax2.set_xlim(-0.8, 14.6)
    ax2.set_ylim(-0.65, len(fams) - 0.35)
    ax2.set_xticks([0, 5, 10])
    ax2.set_xlabel("\u0394pa5 (pp)", labelpad=2)
    fc.clean(ax2)
    ax2.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=fc.GREY, ms=3.4, label="raw"),
                        Line2D([], [], marker="o", ls="", color=fc.BLUE, ms=3.4, label="DiD")],
               loc="lower right", frameon=False, handletextpad=0.25, borderaxespad=0.1,
               labelspacing=0.25, fontsize=8)
    ax2.set_title("5\u2032 1 kb", fontsize=8, pad=3)
    fc.panel(ax2, "d", dx=-0.33, dy=1.10)

    # ---------------- (e) 中下：逐样本差值（自带宽样本轴） ----------------
    ax3 = ax_mm(fig, 87.0, 16.0, 33.0, 46.0)
    loo = ps.did_after_loo.dropna()
    ax3.axvspan(loo.min(), loo.max(), color="#E6E6E6", zorder=0)
    for i, smp in enumerate(order):
        r = ps[ps["sample"] == smp].iloc[0]
        ax3.plot([0, r["diff"]], [i, i], color=fc.LGREY, lw=0.45, zorder=1)
        ax3.plot([r["diff"]], [i], marker="o", ms=3.0, color=fc.colour(smp), zorder=3)
    ax3.axvline(0, color=fc.GREY, lw=0.5)
    ax3.set_yticks(np.arange(ny))
    ax3.set_yticklabels([fc.disp(x) for x in order], fontsize=8)
    ax3.tick_params(axis="y", length=0, pad=2)
    ax3.set_xlim(-1.9, 19.5)
    ax3.set_ylim(ny - 0.45, -0.55)
    ax3.set_xticks([0, 5, 10, 15])
    ax3.set_xlabel("\u0394 vs L1PA4\u20137\n(pp)", labelpad=2)
    fc.clean(ax3)
    ax3.text(loo.mean(), 7.6, f"LOO {loo.min():.1f}\u2013{loo.max():.1f}",
             fontsize=8, ha="center", va="center", color=fc.GREY, rotation=90)
    ax3.set_title("difference", fontsize=8, pad=3)
    fc.panel(ax3, "e", dx=-0.45, dy=1.10)

    # ---------------- (f) 右下：400 bp 窗森林 ----------------
    axc = ax_mm(fig, 139.0, 16.0, 35.0, 46.0)
    famc = ["L1PA2", "L1HS", "LTR7", "LTR5_Hs", "LTR7Y", "L1PA3", "AluY", "AluS",
            "LTR7B", "SVA_D", "SVA_E", "SVA_F"]
    sub = gc[gc.window == "start400"].set_index("family")
    ypos = np.arange(len(famc))[::-1]
    for y, f in zip(ypos, famc):
        r = sub.loc[f]
        axc.errorbar(r.raw_diff, y + 0.16, xerr=[[r.raw_diff - r.raw_lo], [r.raw_hi - r.raw_diff]],
                     fmt="o", ms=2.2, mfc="white", mec=fc.GREY, ecolor=fc.GREY,
                     elinewidth=0.55, capsize=0.9, zorder=3)
        if np.isfinite(r.did):
            axc.errorbar(r.did, y - 0.16, xerr=[[r.did - r.did_lo], [r.did_hi - r.did]],
                         fmt="o", ms=2.6, color=fc.BLUE, ecolor=fc.BLUE, elinewidth=0.55,
                         capsize=0.9, zorder=3)
    axc.axvline(0, color=fc.LGREY, lw=0.5, zorder=1)
    axc.set_yticks(ypos)
    axc.set_yticklabels(famc, fontsize=8)
    axc.set_ylim(-0.7, len(famc) - 0.3)
    axc.set_xlim(-1.4, 13.2)
    axc.set_xticks([0, 4, 8])
    axc.set_xlabel("\u0394pa5 (pp)", labelpad=1)
    fc.clean(axc)
    axc.set_title("400 bp", fontsize=8, pad=3)
    fc.panel(axc, "f", dx=-0.44, dy=1.10)

    # ---------------- 数据与数字 ----------------
    os.makedirs(fc.FIGDIR, exist_ok=True)
    c5.to_csv(os.path.join(fc.FIGDIR, "fig1_data_c5_ranks.tsv"), sep="\t", index=False)
    gc[gc.window.isin(["5p1k", "start400"])].to_csv(
        os.path.join(fc.FIGDIR, "fig1_data_forest.tsv"), sep="\t", index=False)
    ps.to_csv(os.path.join(fc.FIGDIR, "fig1_data_per_sample.tsv"), sep="\t", index=False)

    fc.dump_numbers("fig1", {
        "c5_thresholds": [float(v) for v in c5.threshold],
        "top4_pluripotent": int(c5.top4_is_group.sum()), "n_thresholds": int(len(c5)),
        "intruders": {f"{r.threshold:.2f}": str(r.intruder) for _, r in c5.dropna(subset=["intruder"]).iterrows()},
        "gm27730_ranks": [int(v) for v in c5.rank_GM27730],
        "perm_p_range": [float(c5.perm_p.min()), float(c5.perm_p.max())],
        "pa5_perm_p": float(pa5p.perm_p), "pa5_n_comb": int(pa5p.n_comb),
        "pa5_rank_9": {k: i + 1 for i, k in enumerate(med9.index)},
        "pa5_median_9": {k: round(float(v), 2) for k, v in med9.items()},
        "group_pluri_L1PA2": float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].group_pluri.iloc[0]),
        "group_other_L1PA2": float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].group_other.iloc[0]),
        "raw_L1PA2": [float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].raw_diff.iloc[0]),
                      float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].raw_lo.iloc[0]),
                      float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].raw_hi.iloc[0])],
        "did_L1PA2": [float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].did.iloc[0]),
                      float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].did_lo.iloc[0]),
                      float(gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].did_hi.iloc[0])],
        "per_sample_diff": {r["sample"]: round(float(r["diff"]), 2) for _, r in ps.iterrows()},
        "loo_did": [float(loo.min()), float(loo.max())],
        "start400": {f: [float(sub.loc[f].raw_diff), float(sub.loc[f].did)] for f in famc},
    })

    fc.save(fig, "Fig1", 174)


if __name__ == "__main__":
    main()
