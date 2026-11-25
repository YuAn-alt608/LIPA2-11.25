#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
343 — Fig. 3（174 × 128 mm，5 面板）

(a) 39 个样本的 pa5 中位数（对数）对 r_full，四条子集 ρ
(b) Table 1 的 13 个样本：ρ_obs（空心）与 ρ_corr（实心）+ 噪声带
(c) pa5 与 hiPS6 RNA 的偏相关：实心 = AMB 剔除，空心 = AMB 当未表达
(d) 基因间 403 个拷贝的秩-秩散点，EXP 高亮
(e) H1 Fiber-seq 5mC 与各样本 pa5 的相关（6 点）

(a)(b) 原为 Fig. 2c/2d，因 41 行样本轴在 8 pt 下已占满 Fig. 2，移到这里。

注意：t9m_rna_by_copy.tsv.gz 里的 pa5_* 列是 stride-24 修复之前的旧值，
只能取它的 RNA 列，pa5 一律按 unit_id 从 t9b_pa5_L1PA2_41_fixed 接回。
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

W, H = 174.0, 128.0
SAMPLES = ["H1", "H9", "WTC-11", "GM27730", "HepG2", "Hap1", "K562", "HG002"]


def ax_mm(fig, x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


def main():
    print("Fig. 3")
    fc.spot_check("fig3")

    amb = fc.load("t9m/t9m_partial_rho.tsv")
    allc = fc.load("t13/t13_partial_rho_allcopies.tsv")
    coup = fc.load("t9d/t9d_coupling.tsv")
    rel = fc.rel41()
    rel["r_full"] = pd.to_numeric(rel["r_full"], errors="coerce")
    rel["med_pa5"] = pd.to_numeric(rel["med_pa5"], errors="coerce")

    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))

    # ---------------- (a) 水平 vs 可靠性（原 Fig. 2c） ----------------
    axa = ax_mm(fig, 15.0, 74.0, 60.0, 46.0)
    v = rel[rel.phased.astype(bool)]
    for _, r in v.iterrows():
        c = fc.colour(r["sample"])
        if r.r_s_le0:
            axa.plot([r.med_pa5], [0.0], marker="^", ms=4.0, mfc="white", mec=c, mew=0.8)
        else:
            axa.plot([r.med_pa5], [r.r_full], marker="o", ms=3.4, color=c, mew=0.5)
    axa.set_xscale("log")
    axa.set_xlim(1.0, 22)
    axa.set_ylim(-0.04, 1.0)
    axa.set_xticks([1, 2, 5, 10, 20])
    axa.set_xticklabels(["1", "2", "5", "10", "20"])
    axa.set_yticks([0, 0.5, 1.0])
    axa.set_xlabel("L1PA2 pa5 median (%), log", labelpad=2)
    axa.set_ylabel("r_full", labelpad=2)
    fc.clean(axa)
    pr = fc.load("t15/t15_pa5_vs_rfull.tsv").set_index("subset")
    txt = (f"all 39: ρ = +{pr.loc['all assessed','rho']:.2f} "
           r"($p = 5\times10^{-7}$)" "\n"
           f"non-pluripotent 35: +{pr.loc['non-pluripotent','rho']:.2f}\n"
           f"r_s_nz > 0 (25): +{pr.loc['uncensored r_s_nz>0','rho']:.2f}\n"
           f"LCL 25: +{pr.loc['LCL','rho']:.2f} ($p = 0.33$)")
    axa.text(0.97, 0.04, txt, transform=axa.transAxes, fontsize=8, ha="right", va="bottom")
    fc.style_sample_legend(axa, loc="upper left")
    fc.panel(axa, "a", dx=-0.20, dy=1.015)

    # ---------------- (b) 与 H1 的一致性（原 Fig. 2d） ----------------
    axb = ax_mm(fig, 100.0, 74.0, 68.0, 46.0)
    rho = fc.load("t15/t15_rho_vs_h1_41.tsv", keep_default_na=False, na_values=[""])
    band = fc.load("t10_sim_band.tsv").set_index("sample")
    keep = ["H9", "WTC-11", "GM27730", "Hek293T", "Caco-2", "HepG2", "K562", "MCF-7",
            "Jurkat", "A549", "Panc1", "HG002", "THP-1"]
    rho = rho[rho["sample"].isin(keep)].copy()
    rho["readable"] = rho["readable"].fillna("no")
    rho = rho.sort_values("rho_corr", ascending=True)
    yy = np.arange(len(rho))
    for i, (_, r) in enumerate(rho.iterrows()):
        s = r["sample"]
        c = fc.colour(s)
        if s in band.index:
            b = band.loc[s]
            axb.add_patch(Rectangle((b.band_lo, i - 0.32), b.band_hi - b.band_lo, 0.64,
                                    facecolor="#E6E6E6", edgecolor="none", zorder=0))
        # 不可读的四个样本整行都用灰色（HG002 是 LCL，橙色会与图注"grey"冲突）
        c_obs = c if str(r["readable"]).lower() != "no" else fc.GREY
        axb.plot([r.rho_obs], [i - 0.17], marker="o", ms=3.2, mfc="white", mec=c_obs, mew=0.8, zorder=3)
        if np.isfinite(float(r.rho_corr)):
            if str(r["readable"]).lower() == "no":
                # 不可读：灰色空心 + 灰色虚线 CI，与图注"not readable"一致
                axb.plot([r.lo, r.hi], [i + 0.17, i + 0.17], color=fc.LGREY, lw=0.7,
                         ls=(0, (2, 1.2)), zorder=2)
                axb.plot([r.rho_corr], [i + 0.17], marker="o", ms=3.6, mfc="white",
                         mec=fc.GREY, mew=0.9, zorder=3)
            else:
                axb.plot([r.lo, r.hi], [i + 0.17, i + 0.17], color=c, lw=0.7, zorder=2)
                axb.plot([r.rho_corr], [i + 0.17], marker="o", ms=3.4, color=c, zorder=3)
    axb.set_yticks(yy)
    axb.set_yticklabels(rho["sample"], fontsize=8)
    axb.tick_params(axis="y", length=0, pad=2)
    axb.set_xlim(0, 1.02)
    axb.set_ylim(-0.7, len(rho) - 0.3)
    axb.set_xticks([0, 0.5, 1.0])
    axb.set_xlabel("agreement with H1", labelpad=2)
    fc.clean(axb)
    hd = [Line2D([], [], marker="o", ls="", mfc="white", mec=fc.GREY, mew=0.8, ms=3.4, label="ρ_obs"),
          Line2D([], [], marker="o", ls="", color=fc.GREY, ms=3.4, label="ρ_corr (95% CI)"),
          Line2D([], [], marker="o", ls="", mfc="white", mec=fc.GREY, mew=0.9, ms=3.6,
                 label="ρ_corr, not readable"),
          Rectangle((0, 0), 1, 1, facecolor="#E6E6E6", edgecolor="none", label="noise band")]
    axb.legend(handles=hd, loc="lower right", frameon=False, handletextpad=0.3,
               labelspacing=0.25, fontsize=8, borderaxespad=0.1)
    h9 = rho[rho["sample"] == "H9"].iloc[0]
    b9 = band.loc["H9"]
    axb.annotate(f"+{h9.rho_obs - b9.band_hi:.3f}", (h9.rho_obs, len(rho) - 1 + 0.30),
                 fontsize=8, color=fc.GREY, ha="center", va="bottom")
    fc.panel(axb, "b", dx=-0.12, dy=1.015)

    # ---------------- (c) 偏相关 ----------------
    axc = ax_mm(fig, 15.0, 16.0, 44.0, 44.0)
    ys = np.arange(len(SAMPLES))[::-1]
    C_ALL, C_INT = "#333333", "#009E73"
    for y, s in zip(ys, SAMPLES):
        for scope, colour, dy in (("ALL", C_ALL, 0.26), ("INTERGENIC", C_INT, 0.09)):
            r = amb[(amb.scope == scope) & (amb.pred == s)]
            if len(r):
                r = r.iloc[0]
                axc.plot([r.lo, r.hi], [y + dy, y + dy], color=colour, lw=0.6, zorder=2)
                axc.plot([r.rho_partial], [y + dy], marker="o", ms=3.2, color=colour, zorder=3)
            q = allc[(allc.scope == scope) & (allc.pred == s)]
            if len(q):
                q = q.iloc[0]
                axc.plot([q.lo, q.hi], [y - dy, y - dy], color=colour, lw=0.6, zorder=2,
                         linestyle=(0, (2, 1.2)))
                axc.plot([q.rho_partial], [y - dy], marker="o", ms=3.2, mfc="white",
                         mec=colour, mew=0.8, zorder=3)
    axc.axvline(0, color=fc.LGREY, lw=0.5)
    axc.set_yticks(ys)
    axc.set_yticklabels(SAMPLES, fontsize=8)
    axc.tick_params(axis="y", length=0, pad=2)
    axc.set_xlim(-0.16, 0.68)
    axc.set_ylim(-0.7, len(SAMPLES) + 1.5)   # 顶部留 1.5 行放图例，避免压住 H1/H9
    axc.set_xticks([0, 0.2, 0.4, 0.6])
    axc.set_xlabel("partial ρ (pa5 vs RNA)", labelpad=2)
    fc.clean(axc)
    hc = [Line2D([], [], marker="o", ls="", color=C_ALL, ms=3.2, label="all"),
          Line2D([], [], marker="o", ls="", color=C_INT, ms=3.2, label="intergenic"),
          Line2D([], [], marker="o", ls="", color=fc.GREY, ms=3.2, label="AMB excl."),
          Line2D([], [], marker="o", ls="", mfc="white", mec=fc.GREY, mew=0.8, ms=3.2,
                 label="AMB unexpr.")]
    axc.legend(handles=hc, loc="upper center", frameon=False, ncol=2, handletextpad=0.3,
               labelspacing=0.25, columnspacing=1.2, fontsize=8, borderaxespad=0.0)
    fc.panel(axc, "c", dx=-0.20, dy=1.02)

    # ---------------- (d) RNA 秩-秩散点 ----------------
    axd = ax_mm(fig, 84.0, 16.0, 42.0, 44.0)
    rna = fc.load("t9m/t9m_rna_by_copy.tsv.gz")
    fixed = fc.load("t9b/t9b_pa5_L1PA2_41_fixed.tsv.gz").set_index("unit_id")
    rna["pa5_H1_fixed"] = fixed["H1"].reindex(rna.unit_id).values
    rna["rna"] = np.nanmedian(rna[["senseb_r1", "senseb_r2", "senseb_r3"]].to_numpy(float), axis=1)
    inter = rna[(rna.cat == "intergenic") & rna.label_body.isin(["EXP", "SIL"])].dropna(
        subset=["pa5_H1_fixed", "rna"])
    xr = inter.pa5_H1_fixed.rank(method="average")
    yr = inter.rna.rank(method="average")
    exp = inter.label_body == "EXP"
    axd.plot(xr[~exp], yr[~exp], ls="", marker="o", ms=1.2, mfc=fc.GREY, mec="none", alpha=0.45)
    axd.plot(xr[exp], yr[exp], ls="", marker="o", ms=2.4, color=fc.RED, mec="white", mew=0.3)
    axd.set_xlabel("H1 pa5 rank", labelpad=2)
    axd.set_ylabel("hiPS6 RNA rank", labelpad=2)
    axd.set_xlim(0, len(inter) + 5)
    axd.set_ylim(0, len(inter) + 5)
    axd.set_xticks([0, 200, 400])
    axd.set_yticks([0, 200, 400])
    fc.clean(axd)
    from scipy.stats import spearmanr
    rho_b = float(spearmanr(inter.pa5_H1_fixed, inter.rna).statistic)
    axd.text(0.03, 0.97, f"ρ = {rho_b:.2f}\nn = {len(inter)}\nEXP = {int(exp.sum())}",
             transform=axd.transAxes, fontsize=8, va="top")
    axd.legend(handles=[Line2D([], [], marker="o", ls="", color=fc.RED, ms=5, label="EXP"),
                        Line2D([], [], marker="o", ls="", color=fc.GREY, ms=5, label="SIL")],
               loc="lower right", frameon=False, handletextpad=0.25, labelspacing=0.25,
               fontsize=8, borderaxespad=0.1)
    fc.panel(axd, "d", dx=-0.20, dy=1.02)

    # ---------------- (e) 甲基化（横向，样本名在纵轴） ----------------
    axe = ax_mm(fig, 142.0, 16.0, 32.0, 44.0)
    lab = [f"{r.pa5_of}" for _, r in coup.iterrows()]
    ye = np.arange(len(coup))[::-1]
    for yv, (_, r) in zip(ye, coup.iterrows()):
        same = bool(r.same_cells)
        axe.plot([r.lo, r.hi], [yv, yv], color=fc.GREY, lw=0.7, zorder=2)
        # 菱形=同一批 reads（H1），空心圆=5mC 取自 H1、pa5 取自另一样本；颜色按细胞类型
        col_e = fc.BLUE if same else fc.colour(r.pa5_of)
        axe.plot([r.rho], [yv], marker="D" if same else "o", ms=4.2,
                 color=col_e, mfc=col_e if same else "white",
                 mec=col_e, mew=0.8, zorder=3)
    axe.axvline(0, color=fc.LGREY, lw=0.5)
    axe.set_yticks(ye)
    axe.set_yticklabels(lab, fontsize=8)
    axe.tick_params(axis="y", length=0, pad=2)
    axe.set_xlim(-1.05, 0.10)
    axe.set_xticks([-1.0, -0.5, 0])
    axe.set_ylim(-0.65, len(coup) - 0.35)
    axe.set_xlabel("ρ (5mC vs pa5)", labelpad=2)
    fc.clean(axe)
    axe.annotate("same reads", (coup.rho.iloc[0], ye[0]), xytext=(6, 8),
                 textcoords="offset points", fontsize=8, color=fc.BLUE,
                 va="bottom", ha="left")
    fc.panel(axe, "e", dx=-0.42, dy=1.02)

    # ---------------- 输出 ----------------
    amb.to_csv(os.path.join(fc.FIGDIR, "fig3_data_partial_rho_amb_excluded.tsv"), sep="\t", index=False)
    allc.to_csv(os.path.join(fc.FIGDIR, "fig3_data_partial_rho_amb_negative.tsv"), sep="\t", index=False)
    coup.to_csv(os.path.join(fc.FIGDIR, "fig3_data_meth.tsv"), sep="\t", index=False)
    rho.to_csv(os.path.join(fc.FIGDIR, "fig3_data_rho_h1.tsv"), sep="\t", index=False)
    inter[["unit_id", "cat", "label_body", "pa5_H1_fixed", "rna"]].to_csv(
        os.path.join(fc.FIGDIR, "fig3_data_rna_scatter.tsv"), sep="\t", index=False)
    fc.dump_numbers("fig3", {
        "n_all_amb_excluded": int(amb[amb.scope == "ALL"].n.iloc[0]),
        "n_inter_amb_excluded": int(amb[amb.scope == "INTERGENIC"].n.iloc[0]),
        "h1": {sc: [float(amb[(amb.scope == sc) & (amb.pred == "H1")].rho_partial.iloc[0]),
                    float(amb[(amb.scope == sc) & (amb.pred == "H1")].lo.iloc[0]),
                    float(amb[(amb.scope == sc) & (amb.pred == "H1")].hi.iloc[0])]
               for sc in ("ALL", "INTERGENIC")},
        "k562": {sc: [float(amb[(amb.scope == sc) & (amb.pred == "K562")].rho_partial.iloc[0]),
                      float(amb[(amb.scope == sc) & (amb.pred == "K562")].lo.iloc[0]),
                      float(amb[(amb.scope == sc) & (amb.pred == "K562")].hi.iloc[0])]
                 for sc in ("ALL", "INTERGENIC")},
        "hg002": {sc: [float(amb[(amb.scope == sc) & (amb.pred == "HG002")].rho_partial.iloc[0]),
                       float(amb[(amb.scope == sc) & (amb.pred == "HG002")].lo.iloc[0]),
                       float(amb[(amb.scope == sc) & (amb.pred == "HG002")].hi.iloc[0])]
                  for sc in ("ALL", "INTERGENIC")},
        "h1_amb_negative": {sc: [float(allc[(allc.scope == sc) & (allc.pred == "H1")].rho_partial.iloc[0]),
                                 float(allc[(allc.scope == sc) & (allc.pred == "H1")].lo.iloc[0]),
                                 float(allc[(allc.scope == sc) & (allc.pred == "H1")].hi.iloc[0])]
                            for sc in ("ALL", "INTERGENIC")},
        "rna_scatter": {"n": int(len(inter)), "exp": int(exp.sum()), "rho": rho_b,
                        "zero_frac": float((inter.rna == 0).mean()),
                        "detected": int((inter.rna > 0).sum())},
        "meth": {r.pa5_of: [float(r.rho), float(r.lo), float(r.hi)] for _, r in coup.iterrows()},
        "rho_pa5_rfull": {k: [float(pr.loc[k, "rho"]), float(pr.loc[k, "p"])] for k in pr.index},
        "table1": {r["sample"]: [float(r.rho_obs), float(r.rho_corr), float(r.lo), float(r.hi)]
                   for _, r in rho.iterrows()},
        "h9_above_band": float(h9.rho_obs - b9.band_hi),
    })
    fc.save(fig, "Fig3", 174)


if __name__ == "__main__":
    main()
