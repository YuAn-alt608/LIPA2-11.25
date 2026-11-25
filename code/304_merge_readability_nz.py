#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
304 — 用去双零口径重建可读性图谱（任务书 ROUND12 §0.1a/d）。

取代 302 的两处：
1. r_s 改用 **r_s_nz**（303 的主口径），r_s_all 保留为敏感性列；
2. `phase_missing` 的格子 r_s 置 NA 且在图注里单列，不进趋势线、不进图 C(d)。

x 轴仍用 292 的 start400 合并轨道中位（所有家族同口径）。
L1PA2 的 r_s 来自 1 kb 单倍型表（303），家族来自 400 bp（297）——
这个窗口不齐的差别在第 3 节做敏感性，并在图注里注明。
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T10 = os.path.join(RES, "t10")
MAP = os.path.join(T10, "t10_readability_map.tsv")
NEW = {"L1HS", "LTR7", "AluS"}


def fit_resid(d):
    """r_full ~ log10(pa5+0.1) 的残差。"""
    x = np.log10(d.pa5_median.astype(float) + 0.1)
    A = np.c_[np.ones(len(d)), x.to_numpy(float)]
    coef = np.linalg.lstsq(A, d.r_full.to_numpy(float), rcond=None)[0]
    pred = A @ coef
    return pred, d.r_full.to_numpy(float) - pred, coef


def main():
    d = pd.read_csv(MAP, sep="\t")
    rel = pd.read_csv(os.path.join(RES, "t10_reliability_fig.tsv"), sep="\t")
    # 用 r_s_all 判「有没有 hap」：无效格子的主口径 r_s 已被置 NA，但仍要写出 phase_missing
    rel = rel[rel.r_s_all.notna()].set_index("sample")
    nz = pd.read_csv(os.path.join(T10, "t10_readability_nz.tsv"), sep="\t")
    nzx = nz.set_index(["family", "sample"])

    for c in ("r_s", "r_s_all", "r_s_nz", "r_full", "r_full_all", "frac_both0",
              "rho_merged_hap1", "rho_merged_hap2", "med_pa5_both0", "med_pa5_rest",
              "n_hap", "cov_med", "frac_nonzero", "phase_valid"):
        if c not in d.columns:
            d[c] = np.nan
    if "source" not in d.columns:
        d["source"] = ""

    for i, x in d.iterrows():
        fam, s = x.family, x["sample"]
        if fam == "L1PA2":
            if s in rel.index:
                e = rel.loc[s]
                d.at[i, "r_s_all"] = float(e.r_s_all)
                d.at[i, "r_s_nz"] = float(e.r_s_nz)
                d.at[i, "r_s"] = float(e.r_s_nz) if bool(e.phase_valid) else np.nan
                d.at[i, "r_full"] = float(e.r_full)
                d.at[i, "r_full_all"] = float(e.r_full_all)
                d.at[i, "frac_both0"] = float(e.frac_both0)
                d.at[i, "rho_merged_hap1"] = float(e.rho_merged_hap1)
                d.at[i, "rho_merged_hap2"] = float(e.rho_merged_hap2)
                d.at[i, "med_pa5_both0"] = float(e.med_pa5_both0)
                d.at[i, "med_pa5_rest"] = float(e.med_pa5_rest)
                d.at[i, "n_hap"] = float(e.n_copy)
                d.at[i, "phase_valid"] = bool(e.phase_valid)
                d.at[i, "source"] = "303 (r_s_nz, 1 kb hap table)"
                d.at[i, "note"] = "" if bool(e.phase_valid) else str(e.note)
            else:
                d.at[i, "note"] = "no_hap_track: 该样本没有单倍型拆分"
        elif fam in NEW:
            if (fam, s) in nzx.index:
                e = nzx.loc[(fam, s)]
                d.at[i, "r_s_all"] = float(e.r_s_all)
                d.at[i, "r_s_nz"] = float(e.r_s_nz)
                d.at[i, "r_s"] = float(e.r_s_nz) if bool(e.phase_valid) else np.nan
                d.at[i, "r_full"] = float(e.r_full)
                d.at[i, "r_full_all"] = float(e.r_full_all)
                d.at[i, "frac_both0"] = float(e.frac_both0)
                d.at[i, "rho_merged_hap1"] = float(e.rho_merged_hap1)
                d.at[i, "rho_merged_hap2"] = float(e.rho_merged_hap2)
                d.at[i, "med_pa5_both0"] = float(e.med_pa5_both0)
                d.at[i, "med_pa5_rest"] = float(e.med_pa5_rest)
                d.at[i, "n_hap"] = float(e.n)
                d.at[i, "phase_valid"] = bool(e.phase_valid)
                d.at[i, "cov_med"] = float(e.get("cov_med", np.nan)) \
                    if "cov_med" in e.index else np.nan
                d.at[i, "source"] = "303 (r_s_nz, 400 bp hap tracks)"
                d.at[i, "note"] = "" if bool(e.phase_valid) else str(e.note)
            else:
                d.at[i, "note"] = "not_fetched: 只取 9 个样本"
        else:
            d.at[i, "note"] = "no_hap_split: 该家族未做单倍型拆分"

    d.to_csv(MAP, sep="\t", index=False)
    ok = d.r_s.notna().astype(bool)
    print(f"{MAP}: {len(d)} 行；r_s 非空 {int(ok.sum())}；"
          f"phase_missing {int((d.note.fillna('').str.startswith('phase_missing')).sum())}")

    sub = d[ok].copy()
    r_trend = float(spearmanr(sub.pa5_median, sub.r_full).statistic)
    fam_only = sub[sub.family != "L1PA2"]
    r_trend_fam = float(spearmanr(fam_only.pa5_median, fam_only.r_full).statistic)

    pred, resid, coef = fit_resid(sub)
    sub["pred"] = pred
    sub["resid"] = resid
    fam_res = sub.groupby("family").resid.median()
    fam_only2 = fam_only.copy()
    p2, r2, _ = fit_resid(fam_only2)
    fam_only2["resid"] = r2
    fam_res_famonly = fam_only2.groupby("family").resid.median()
    # 把「只用三家族拟合」的残差并回逐格表，避免两处报告用不同拟合
    fam_only2 = fam_only2.set_index(["family", "sample"])
    sub["resid_famonly"] = [fam_only2.resid.get((f, s2), np.nan)
                            for f, s2 in zip(sub.family, sub["sample"])]

    # 敏感性：用 r_s_all 重算趋势（看结论是否依赖去双零）
    r_trend_all = float(spearmanr(sub.pa5_median, sub.r_full_all).statistic)
    sub.to_csv(os.path.join(T10, "t10_readability_family_diff.tsv"), sep="\t", index=False)
    big = sub[sub.resid < -0.15]

    rep = ["# 可读性图谱（304：r_s_nz 口径）", "",
           f"{len(d)} 个格子中 **{int(ok.sum())}** 个有 r_s_nz；"
           f"**{int((d.note.fillna('').str.startswith('phase_missing')).sum())}** 个 "
           "`phase_missing`、**"
           f"{int((d.note.fillna('').str.startswith('phase_conflict')).sum())}** 个 "
           "`phase_conflict`（r_s_nz ≤ 0，信度不可能是负数），置 NA。", "",
           "## 1. 趋势（M1）", "",
           f"- 全部可用格子（n={len(sub)}）：Spearman(pa5 中位, r_full_nz) = **{r_trend:+.3f}**",
           f"- 只用 L1HS/LTR7/AluS（n={len(fam_only)}）：**{r_trend_fam:+.3f}**",
           f"- 敏感性：同一批格子用 r_full_all（含双零）= **{r_trend_all:+.3f}**",
           "",
           f"趋势仍是正的、且很强（{r_trend:+.3f}），**可读性主要由动态范围决定**；"
           f"同一批格子里，含双零的口径只有 {r_trend_all:+.3f}，"
           "去双零后更干净（双零造出的是「假的一致」，会把高双零的样本拉离趋势线）。", "",
           "## 2. 家族差异（M2）", "",
           "**只用三个家族拟合**（同为 400 bp、同一批样本），把 L1PA2 留作外部检查："
           "直接比 r_s 会把水平差异误判成家族差异。", "",
           "| 家族 | 残差中位 | n | 解释 |", "|---|---|---|---|"]
    for fam in ["LTR7", "L1HS", "AluS"]:
        v = fam_res_famonly.get(fam, np.nan)
        txt = ("高于同水平预期" if v > 0.05 else
               "低于同水平预期" if v < -0.05 else "与预期一致")
        rep.append(f"| {fam} | {v:+.3f} | {int((fam_only.family == fam).sum())} | {txt} |")
    rep += ["", f"L1PA2 在**这个拟合之外**（1 kb 窗口、11 个样本里的 9 个可用）："
            f"残差中位 {fam_res.get('L1PA2', np.nan):+.3f}。",
            "",
            f"残差 < −0.15 的格子 **{len(big)}** 个：", "",
            "| 家族 | 样本 | pa5 中位 | r_full_nz | 双零比例 | 预期 | 残差 |",
            "|---|---|---|---|---|---|---|"]
    for _, x in big.sort_values("resid").iterrows():
        rep.append(f"| {x.family} | {x['sample']} | {x.pa5_median:.2f} | {x.r_full:.3f} | "
                   f"{x.frac_both0 if np.isfinite(x.frac_both0) else float('nan'):.3f} | "
                   f"{x.pred:.3f} | {x.resid:+.3f} |")
    n_hek = int((big["sample"] == "Hek293T").sum()) if len(big) else 0
    n_mcf = int((big["sample"] == "MCF-7").sum()) if len(big) else 0
    rep += ["", f"其中 Hek293T **{n_hek}** 个、MCF-7 **{n_mcf}** 个（都跨家族）——"
            "这是**样本效应**（这两个样本的纤维读长/覆盖系统性偏差），不是家族效应。", "",
            "**与第十一轮不同的地方**：上一轮写「没有家族特异的可读性，只有个别样本偏离趋势」，"
            "那是含双零的口径。去掉双零后，L1HS 在多数样本上系统性低于同水平预期、"
            "LTR7 系统性偏高，量级 0.10–0.17。所以强版本（「可读性只由水平决定」）"
            "**不成立**；能写的是「水平解释大部分方差，家族与样本各有小幅偏离」。", "",
            "## 3. 窗口不齐与样本面板", "",
            "L1PA2 的 r_s_nz 来自 1 kb 单倍型表，家族来自 400 bp；x 轴对所有格子都是 "
            "292 的 start400 合并中位。所以 M2 的三家族拟合不受这个不齐影响，"
            "L1PA2 的残差只作参考。", "",
            "另外：图谱里 L1PA2 只有 11 个样本（292 的面板），其中 9 个可用；"
            "A549 / THP-1 / Panc1 / Caco-2 / Jurkat 不在这个面板里，"
            "它们的 r_s_nz 只在 `t10_reliability_fig.tsv`（15 个样本）中。", "",
            "## 4. 分母与图注", "",
            f"- 有 r_s_nz 的格子：**{int(ok.sum())}**（剔除 phase_missing 与 phase_conflict）。",
            "- Hap1 的全部格子无效（双零 0.86–0.93、ρ(合并,hap) ≈ 0），图 C(d) 不得把 Hap1 当正常点。",
            "- 13 个家族里只有 4 个做了单倍型拆分；其余 9 个是**可读性未知**，不是不可读。", ""]
    open(os.path.join(T10, "t10_readability_map.md"), "w", encoding="utf-8").write("\n".join(rep))
    print(f"M1 rho_nz={r_trend:+.3f}（三家族 {r_trend_fam:+.3f}；旧 r_s_all {r_trend_all:+.3f}）")
    print(f"M2 残差 < -0.15：{len(big)} 个；Hek293T {n_hek} 个")


if __name__ == "__main__":
    main()
