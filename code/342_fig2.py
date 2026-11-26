#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
342 — Fig. 2（174 × 185 mm，41 个样本的可靠性双联）

(a) 41 样本的三个 r_s：r_s_all（含双零）、r_s_nz（去双零）、r_s_cov（覆盖感知，
    仅 25 个可评估且有覆盖轨道的样本）
(b) 39 个通过定相样本的 r_full + 95% 区间，0.667 阈值；r_s_nz<=0 画在 0 的空心三角

8 pt 字下 41 行至少需要 ~127 mm 面板高，所以 (c)(d) 已移到 Fig. 3，
本图只保留样本轴，图例放到画布底部。
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

W, H = 174.0, 185.0
THR = 0.667


def ax_mm(fig, x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


def main():
    print("Fig. 2")
    fc.spot_check("fig2")

    rel = fc.rel41()
    for c in ("r_full", "r_s_all", "r_s_nz", "r_full_lo", "r_full_hi", "med_pa5"):
        rel[c] = pd.to_numeric(rel[c], errors="coerce")

    cov = pd.read_csv(os.path.join(fc.RES, "t15", "t15_table_s4_cov.tsv"), sep="\t")
    cov = cov.set_index("sample")[["n_cov", "frac_both0_covered", "r_s_cov", "r_full_cov"]]
    rel = rel.join(cov, on="sample")
    for c in ("r_s_cov", "r_full_cov", "frac_both0_covered"):
        rel[c] = pd.to_numeric(rel[c], errors="coerce")

    # r28：主估计量换成覆盖感知版。把 389 的 r_full_cov 与区间顶到 r_full 系列上，
    # 旧的双零口径留成 r_full_nz（图 2b 不再画它，数字仍在 fig2_numbers.json 里可查）。
    new = fc.load("t15/t15_reliability_cov_41.tsv").set_index("sample")
    rel = rel.join(new[["r_full_cov_lo", "r_full_cov_hi", "readable"]].rename(
        columns={"readable": "readable_cov"}), on="sample")
    rel["r_full_nz"] = rel["r_full"]
    rel["r_full_lo_nz"] = rel["r_full_lo"]
    rel["r_full_hi_nz"] = rel["r_full_hi"]
    rel["r_full"] = rel["r_full_cov"]
    rel["r_full_lo"] = rel["r_full_cov_lo"]
    rel["r_full_hi"] = rel["r_full_cov_hi"]
    rel["readable"] = rel["readable_cov"]

    order = fc.sample_order_by_rfull(rel)
    rel = rel.set_index("sample").loc[order].reset_index()
    y = np.arange(len(rel))[::-1]
    rel["y"] = y
    ny = len(rel)

    fig = plt.figure(figsize=(W * fc.MM, H * fc.MM))

    # ---------------- (a) 哑铃 ----------------
    axa = ax_mm(fig, 40.0, 36.0, 52.0, 140.0)
    for _, r in rel.iterrows():
        if not bool(r.phased):
            axa.plot([0.02], [r.y], marker="x", ms=4.0, color=fc.LGREY, mew=0.9)
            continue
        axa.plot([r.r_s_all, r.r_s_nz], [r.y, r.y], color=fc.LGREY, lw=0.6, zorder=1)
        axa.plot([r.r_s_all], [r.y], marker="o", ms=3.4, mfc="white",
                 mec=fc.colour(r["sample"]), mew=0.8, zorder=3)
        axa.plot([r.r_s_nz], [r.y], marker="o", ms=3.4, color=fc.colour(r["sample"]), zorder=3)
        if np.isfinite(r.r_s_cov):
            axa.plot([r.r_s_cov], [r.y], marker="s", ms=2.9, mfc="#222222", mec="none", zorder=4)
    axa.axvline(0, color=fc.LGREY, lw=0.5)
    axa.set_xlim(-0.45, 1.02)
    axa.set_ylim(-0.8, ny - 0.2)
    axa.set_yticks(y)
    axa.set_yticklabels(rel["sample"], fontsize=8)
    axa.tick_params(axis="y", length=0, pad=2)
    axa.set_xticks([0, 0.5, 1.0])
    axa.set_xlabel("split-half correlation r_s", labelpad=2)
    fc.clean(axa)
    fc.panel(axa, "a", dx=-0.66, dy=1.002)

    # ---------------- (b) r_full ----------------
    axb = ax_mm(fig, 104.0, 36.0, 66.0, 140.0)
    for _, r in rel.iterrows():
        if not bool(r.phased):
            axb.plot([0.0], [r.y], marker="x", ms=4.0, color=fc.LGREY, mew=0.9)
            continue
        c = fc.colour(r["sample"])
        style = dict(marker="o", ms=4.0, color=c, zorder=3)
        if r.readable == "borderline":
            style.update(mfc="white", mfcalt=c, fillstyle="left", mew=0.7)
        elif r.readable == "no":
            style.update(mfc="white", mew=0.7)
        if np.isfinite(r.r_full_lo):
            axb.plot([r.r_full_lo, r.r_full_hi], [r.y, r.y], color=c, lw=0.7, zorder=2)
        axb.plot([r.r_full], [r.y], **style)
    axb.axvline(THR, color=fc.GREY, lw=0.7, ls=(0, (3, 2)))
    axb.text(THR + 0.015, ny - 0.4, "0.667", fontsize=8, color=fc.GREY, va="top")
    axb.set_xlim(-0.05, 1.0)
    axb.set_ylim(-0.8, ny - 0.2)
    axb.set_yticks([])
    axb.set_xticks([0, 0.5, 1.0])
    axb.set_xlabel("r_full,cov (Spearman–Brown)", labelpad=2)
    fc.clean(axb, left=False)
    fc.panel(axb, "b", dx=-0.12, dy=1.002)

    # ---------------- 图例（画布底部，两栏） ----------------
    ha = [Line2D([], [], marker="o", ls="", mfc="white", mec=fc.GREY, mew=0.8, ms=3.6,
                 label="with double-zero"),
          Line2D([], [], marker="o", ls="", color=fc.GREY, ms=3.6, label="double-zero removed"),
          Line2D([], [], marker="s", ls="", mfc="#222222", mec="none", ms=3.0,
                 label="coverage-aware"),
          Line2D([], [], marker="x", ls="", color=fc.LGREY, mew=0.9, ms=3.6, label="not assessed")]
    hb = [Line2D([], [], marker="o", ls="", color=fc.GREY, ms=3.6, label="readable"),
          Line2D([], [], marker="o", ls="", mfc="white", mfcalt=fc.GREY, fillstyle="left",
                 mec=fc.GREY, mew=0.7, ms=3.6, label="borderline"),
          Line2D([], [], marker="o", ls="", mfc="white", mec=fc.GREY, mew=0.7, ms=3.6,
                 label="not readable")]

    # 两排图例都用 mode="expand" 铺满同一条左右边线：每项间距相等，
    # 不会因为某一排项数少（3 项）就挤到一侧或折成两行。
    for handles, ncol, ybot in ((ha, 4, 0.058), (hb, 3, 0.016)):
        fig.legend(handles=handles, loc="lower left",
                   bbox_to_anchor=(0.0, ybot, 0.992, 0.001),
                   ncol=ncol, mode="expand", frameon=False, borderpad=0,
                   handletextpad=0.3, handlelength=1.2, columnspacing=1.6,
                   fontsize=8, borderaxespad=0)

    # ---------------- (c)(d) 的数（画在 Fig. 3，此处只登记） ----------------
    pr = fc.load("t15/t15_pa5_vs_rfull_cov.tsv").set_index("subset")
    rho = fc.load("t15/t15_rho_vs_h1_cov_41.tsv", keep_default_na=False, na_values=[""])
    band = fc.load("t10_sim_band_cov.tsv").set_index("sample")
    keep = ["H9", "WTC-11", "GM27730", "Hek293T", "Caco-2", "HepG2", "K562", "MCF-7",
            "Jurkat", "A549", "Panc1", "HG002", "THP-1"]
    rho = rho[rho["sample"].isin(keep)].copy()
    rho["readable"] = rho["readable"].fillna("no")
    rho = rho.sort_values("rho_corr", ascending=True)
    h9 = rho[rho["sample"] == "H9"].iloc[0]
    b9 = band.loc["H9"]

    # ---------------- 输出 ----------------
    rel.to_csv(os.path.join(fc.FIGDIR, "fig2_data_reliability.tsv"), sep="\t", index=False)
    fc.dump_numbers("fig2", {
        "n_total": int(len(rel)), "n_phased": int(rel.phased.sum()),
        "n_readable": int(rel.readable.isin(["yes", "borderline"]).sum()),
        "n_readable_clear": int((rel.readable == "yes").sum()),
        "n_borderline": int((rel.readable == "borderline").sum()),
        "n_not_readable": int((rel.readable == "no").sum()),
        "n_not_assessed": int((rel.readable == "n/a").sum()),
        "n_rs_le0": int(rel.r_s_le0.sum()),
        "cov": (lambda ca: {"n": int(rel.r_s_cov.notna().sum()), "n_assessable": len(ca),
                            "r_s_min": float(ca.r_s_cov.min()), "r_s_max": float(ca.r_s_cov.max()),
                            "n_le0": int((ca.r_s_cov <= 0).sum()),
                            "n_ge_half": int((ca.r_s_cov >= 0.5).sum()),
                            "r_full_max": float(ca.r_full_cov.max()),
                            "frac_both0_cov_min": float(ca.frac_both0_covered.min()),
                            "frac_both0_cov_max": float(ca.frac_both0_covered.max())})(
            rel[rel.r_s_cov.notna() & rel.phased]),
        "hap1": {"r_s_all": float(rel.set_index("sample").loc["Hap1", "r_s_all"]),
                 "r_s_nz": float(rel.set_index("sample").loc["Hap1", "r_s_nz"])},
        "panc1": [float(rel.set_index("sample").loc["Panc1", "r_s_all"]),
                  float(rel.set_index("sample").loc["Panc1", "r_s_nz"])],
        "a549": [float(rel.set_index("sample").loc["A549", "r_s_all"]),
                 float(rel.set_index("sample").loc["A549", "r_s_nz"])],
        "thp1": [float(rel.set_index("sample").loc["THP-1", "r_s_all"]),
                 float(rel.set_index("sample").loc["THP-1", "r_s_nz"])],
        # 下列两组数现在画在 Fig. 3（原 Fig. 2c/2d），仍在此登记以便审计沿用同一入口
        "rho_pa5_rfull": {k: [float(pr.loc[k, "rho"]), float(pr.loc[k, "p"])] for k in pr.index},
        "table1": {r["sample"]: [float(r.rho_obs), float(r.rho_corr), float(r.lo), float(r.hi)]
                   for _, r in rho.iterrows()},
        "h9_above_band": float(h9.rho_obs - b9.band_hi),
    })
    fc.save(fig, "Fig2", 174)


if __name__ == "__main__":
    main()
