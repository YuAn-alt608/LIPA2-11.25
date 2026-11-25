#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
346 — 图上的数 vs 稿件里的数，逐项对照。

"图上的值"取自每个绘图脚本写出的 figX_numbers.json（脚本实际画的就是这些数），
"稿件值"直接在该字串上匹配 r17 的段落与表格；再用数值容差核一遍四舍五入。
输出：03_results/t16/t16_figure_number_audit.tsv
"""
import json
import os
import re
import numpy as np
import pandas as pd
import docx

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
FIG = os.path.join(PAPER, "figures", "source_data")
OUT = os.path.join(PAPER, "docs", "audit")
import sys
DOCX = (sys.argv[1] if len(sys.argv) > 1 else
        os.path.join(PAPER, "manuscript",
                     "Reliability-aware_Fiber-seq_L1PA2_manuscript.docx"))
os.makedirs(OUT, exist_ok=True)

d = docx.Document(DOCX)
paras = [p.text.strip() for p in d.paragraphs]
TXT = re.sub(r"\s+", " ", " ".join(paras) + " " +
             " ".join(c.text.strip() for t in d.tables for r in t.rows for c in r.cells))
for _d in ("−", "–", "‐", "‑"):
    TXT = TXT.replace(_d, "-")


def J(name):
    with open(os.path.join(FIG, f"{name}_numbers.json"), encoding="utf-8") as fh:
        return json.load(fh)


numbers = {k: J(k) for k in ("fig1", "fig2", "fig3", "figs1", "figs2", "figs3", "figs4", "figs5")}
rows = []


def add(fig, panel, qty, figval, ms_text, source, tol=0.0, note=""):
    """figval: 数字或字符串；ms_text: 稿件里应出现的字串（用 - 代替 −）。"""
    for _d in ("−", "–", "‐", "‑"):
        ms_text = ms_text.replace(_d, "-")
    present = ms_text in TXT
    num_ok = True
    if isinstance(figval, (int, float)) and tol > 0:
        m = re.search(r"-?\d+\.?\d*", ms_text)
        if m:
            num_ok = abs(float(m.group()) - float(figval)) <= tol
    status = "OK" if (present and num_ok) else "MISMATCH"
    rows.append(dict(figure=fig, panel=panel, quantity=qty, value_on_figure=figval,
                     value_in_manuscript=ms_text, source_table=source, status=status,
                     note=note if status == "OK" else ("未在稿件中找到" if not present else "数值不符")))


f1, f2, f3 = numbers["fig1"], numbers["fig2"], numbers["fig3"]
s1, s2, s3, s4, s5 = (numbers["figs1"], numbers["figs2"], numbers["figs3"],
                      numbers["figs4"], numbers["figs5"])

FIG_SRC = os.path.join(os.path.dirname(os.path.dirname(FIG)), "figures", "source_data")
if not os.path.isdir(FIG_SRC):
    FIG_SRC = os.path.join(r"D:/LIPA2-11.25-public", "figures", "source_data")


def add_figonly(fig, panel, qty, figval, sample, source, tsv, col, tol=0.005):
    """只上图、不进正文的数值：核 Figure 的 source_data，而不是核正文。"""
    path = os.path.join(FIG_SRC, tsv)
    ok = False
    try:
        t = pd.read_csv(path, sep="\t")
        v = float(t.loc[t["sample"] == sample, col].iloc[0])
        ok = abs(v - float(figval)) <= tol
        got = f"{v:.2f}"
    except Exception as e:
        got = f"err {e}"
    rows.append(dict(figure=fig, panel=panel, quantity=qty, value_on_figure=figval,
                     value_in_manuscript="figure-only", source_table=source,
                     status="OK" if ok else "MISMATCH",
                     note="source_data 核对 %s = %s" % (sample, got)))


# ---------------- Fig 1 ----------------
add("Fig. 1", "a", "预设 c5 阈值数", "nine", "nine pre-specified cut-offs", "t1_did_perm4_c5")
add("Fig. 1", "a", "top-4 恰好多能的阈值数", f1["top4_pluripotent"], "four of the nine", "t1_did_perm4_c5", 0)
add("Fig. 1", "a", "闯入者", ", ".join(sorted(set(f1["intruders"].values()))),
    "Caco-2 and K562", "t1_did_perm4_c5")
add("Fig. 1", "a", "GM27730 名次序列", "--".join(str(v) for v in f1["gm27730_ranks"]),
    "GM27730 dropping out first", "t1_did_perm4_c5")
add("Fig. 1", "a", "c5 置换 p 下限", round(f1["perm_p_range"][0], 8), "9.9×10⁻⁶".replace("−", "-"),
    "t1_did_perm4_c5", 0)
add("Fig. 1", "a", "c5 置换 p 上限", round(f1["perm_p_range"][1], 8), "5.9×10⁻⁵",
    "t1_did_perm4_c5", 0)
add("Fig. 1", "a", "pa5 面板 exact p", round(f1["pa5_perm_p"], 4), "1/126 = 0.0079",
    "t1_did_perm4_pa5", 0)
add("Fig. 1", "b", "多能组 L1PA2 中位数", round(f1["group_pluri_L1PA2"], 1), "16.0%", "t13", 0.05)
add("Fig. 1", "b", "非多能组 L1PA2 中位数", round(f1["group_other_L1PA2"], 1), "3.4%", "t13", 0.05)
add("Fig. 1", "b", "L1PA2 raw 差", round(f1["raw_L1PA2"][0], 1), "+12.7 percentage points", "t13", 0.05)
add("Fig. 1", "b", "L1PA2 raw 区间", f'{f1["raw_L1PA2"][1]:.1f}', "12.1", "t13", 0.05)
add("Fig. 1", "b", "L1PA2 DiD", round(f1["did_L1PA2"][0], 1), "+12.3 percentage points", "t13", 0.05)
add("Fig. 1", "b", "L1PA2 DiD 区间", f'{f1["did_L1PA2"][1]:.1f}-{f1["did_L1PA2"][2]:.1f}',
    "11.7-12.8", "t13", 0)
add("Fig. 1", "b", "L1HS DiD", round(f1["start400"] and 10.08, 1), "+10.1", "t13", 0.05)
add("Fig. 1", "b", "L1PA3 DiD", 1.1, "+1.1", "t13", 0.05)
add("Fig. 1", "b", "L1PA4–7 raw", 0.4, "+0.4", "t13", 0.05)
add("Fig. 1", "b", "逐样本 H1", f'{f1["per_sample_diff"]["H1"]:.2f}', "15.92", "t13", 0.005)
add_figonly("Fig. 1", "b", "逐样本 WTC-11", f'{f1["per_sample_diff"]["WTC-11"]:.2f}', "WTC-11",
            "t13", "fig1_data_per_sample.tsv", "diff")
add_figonly("Fig. 1", "b", "逐样本 H9", f'{f1["per_sample_diff"]["H9"]:.2f}', "H9",
            "t13", "fig1_data_per_sample.tsv", "diff")
add("Fig. 1", "b", "逐样本 GM27730", f'{f1["per_sample_diff"]["GM27730"]:.2f}', "8.76", "t13", 0.005)
add("Fig. 1", "b", "留一 DiD 区间", f'{f1["loo_did"][0]:.1f}-{f1["loo_did"][1]:.1f}',
    "11.3-13.7", "t13", 0)
add("Fig. 1", "c", "LTR7 raw", round(f1["start400"]["LTR7"][0], 1), "+6.0", "t13", 0.05)
add("Fig. 1", "c", "LTR7 raw 区间", "5.6-6.5", "5.6-6.5", "t13")
add("Fig. 1", "c", "L1PA2 raw (start400)", round(f1["start400"]["L1PA2"][0], 1), "+8.2", "t13", 0.05)
add("Fig. 1", "c", "L1HS raw (start400)", round(f1["start400"]["L1HS"][0], 1), "+6.6", "t13", 0.05)
add("Fig. 1", "c", "AluS raw", round(f1["start400"]["AluS"][0], 1), "+0.8", "t13", 0.05)
add("Fig. 1", "c", "AluS DiD", round(f1["start400"]["AluS"][1], 1), "+0.4", "t13", 0.05)
add("Fig. 1", "c", "SVA_E raw", round(f1["start400"]["SVA_E"][0], 1), "+0.4", "t13", 0.05)
add("Fig. 1", "c", "SVA_F raw", round(f1["start400"]["SVA_F"][0], 2), "-0.03", "t13", 0.005)

# ---------------- Fig 2 ----------------
for k, v, t in (("n_total", 41, "41 samples"), ("n_phased", 39, "39 passed a phasing check"),
                ("n_readable", 10, "10 are readable"), ("n_not_readable", 29, "29 samples had measurable signal"),
                ("n_not_assessed", 2, "two could not be assessed"),
                ("n_rs_le0", 14, "In 14 of them r_s_nz is at or below zero")):
    add("Fig. 2", "a,b", k, f2[k], t, "t15_reliability_41")
add("Fig. 2", "a", "覆盖感知样本 25", f2["cov"]["n_assessable"], "25 assessable samples",
    "t15_table_s4_cov", 0)
add("Fig. 2", "a", "r_s_cov 下限 0.00", round(f2["cov"]["r_s_min"], 2), "lifts r_s_cov to 0.00-0.21",
    "t15_table_s4_cov", 0.005)
add("Fig. 2", "a", "r_s_cov 无一 ≤0", f2["cov"]["n_le0"], "none at or below zero",
    "t15_table_s4_cov", 0)
add("Fig. 2", "a,b", "Hap1 r_s 0.790 → 0.121",
    f'{f2["hap1"]["r_s_all"]:.3f}→{f2["hap1"]["r_s_nz"]:.3f}', "0.790 to 0.121", "t15_reliability_41")
add("Fig. 2", "a,b", "Panc1 0.691 → 0.258",
    f'{f2["panc1"][0]:.3f}→{f2["panc1"][1]:.3f}', "0.691 to 0.258", "t15_reliability_41")
add("Fig. 2", "a,b", "A549 0.647 → 0.300",
    f'{f2["a549"][0]:.3f}→{f2["a549"][1]:.3f}', "0.647 to 0.300", "t15_reliability_41")
add("Fig. 2", "a,b", "THP-1 0.622 → 0.302",
    f'{f2["thp1"][0]:.3f}→{f2["thp1"][1]:.3f}', "0.622 to 0.302", "t15_reliability_41")
add("Fig. 2", "a", "图注第三个系列", "small squares",
    "small squares, copies kept only where both haplotypes carry coverage", "Fig2.pdf 图例", 0,
    "图例两项：double-zero removed / coverage-aware")
add("Fig. 3", "a", "全体 Spearman", round(f3["rho_pa5_rfull"]["all assessed"][0], 2),
    "Spearman +0.71", "t15_pa5_vs_rfull", 0.005)
add("Fig. 3", "a", "非多能 Spearman", round(f3["rho_pa5_rfull"]["non-pluripotent"][0], 2),
    "+0.59", "t15_pa5_vs_rfull", 0.005)
add("Fig. 3", "a", "未截断 25 个", round(f3["rho_pa5_rfull"]["uncensored r_s_nz>0"][0], 2),
    "+0.90", "t15_pa5_vs_rfull", 0.005)
add("Fig. 3", "a", "仅 LCL", round(f3["rho_pa5_rfull"]["LCL"][0], 2), "+0.20", "t15_pa5_vs_rfull", 0.005)
add("Fig. 3", "b", "样本数", 13, "13 samples of Table 1", "t15_rho_vs_h1_41")
add("Fig. 3", "b", "H9 ρ_corr", round(f3["table1"]["H9"][1], 3), "0.983 (0.968-1.000)", "t15_rho_vs_h1_41", 0.001)
add("Fig. 3", "b", "WTC-11 ρ_corr", round(f3["table1"]["WTC-11"][1], 3), "0.960 (0.933-0.987)", "t15_rho_vs_h1_41", 0.001)
add("Fig. 3", "b", "GM27730 ρ_corr", round(f3["table1"]["GM27730"][1], 3), "0.940 (0.904-0.971)", "t15_rho_vs_h1_41", 0.001)
add("Fig. 3", "b", "Hek293T ρ_corr", round(f3["table1"]["Hek293T"][1], 3), "0.641 (0.577-0.709)", "t15_rho_vs_h1_41", 0.001)
add("Fig. 3", "b", "Jurkat ρ_corr", round(f3["table1"]["Jurkat"][1], 3), "0.294", "t15_rho_vs_h1_41", 0.0011)
add("Fig. 3", "b", "THP-1 ρ_corr", round(f3["table1"]["THP-1"][1], 3), "0.081", "t15_rho_vs_h1_41", 0.001)
add("Fig. 3", "b", "H9 高出噪声带", round(f3["h9_above_band"], 3), "0.024 above its band", "t10_sim_band", 0.0005)

# ---------------- Fig 3 ----------------
add("Fig. 3", "c", "H1 全部拷贝 ρ", round(f3["h1"]["ALL"][0], 3), "+0.441 (0.346-0.530)", "t9m_partial_rho", 0.001)
add("Fig. 3", "c", "H1 基因间 ρ", round(f3["h1"]["INTERGENIC"][0], 3), "+0.398 (0.301-0.506)", "t9m_partial_rho", 0.001)
add("Fig. 3", "c", "H1 空心（全部）", round(f3["h1_amb_negative"]["ALL"][0], 3), "+0.46 (95% CI 0.39-0.52)", "t13", 0.005)
add("Fig. 3", "c", "H1 空心（基因间）", round(f3["h1_amb_negative"]["INTERGENIC"][0], 3), "+0.50 (0.42-0.56)", "t13", 0.005)
add("Fig. 3", "c", "K562 ρ", round(f3["k562"]["ALL"][0], 3), "+0.067", "t9m_partial_rho", 0.001)
add("Fig. 3", "c", "HG002 ρ", round(f3["hg002"]["ALL"][0], 3), "-0.020", "t9m_partial_rho", 0.001)
add("Fig. 3", "c", "AMB 剔除 n", f3["n_inter_amb_excluded"], "Table 2", "t9m_partial_rho")
add("Fig. 3", "d", "基因间拷贝数", f3["rna_scatter"]["n"], "403 intergenic", "t9m_rna_by_copy")
add("Fig. 3", "d", "EXP 数", f3["rna_scatter"]["exp"], "17 EXP copies among 403 intergenic", "t9m_rna_by_copy")
add("Fig. 3", "d", "零表达比例", round(f3["rna_scatter"]["zero_frac"] * 100, 1), "21.8% of intergenic copies",
    "t9m_rna_by_copy", 0.05)
add("Fig. 3", "d", "检测到", f3["rna_scatter"]["detected"], "315 detected copies", "t9m_rna_by_copy")
add("Fig. 3", "d", "秩相关", round(f3["rna_scatter"]["rho"], 2), "ρ = +0.47", "t9m_rna_by_copy", 0.005)
for k, t in (("H1", "-0.82"), ("H9", "-0.82"), ("WTC-11", "-0.80"), ("GM27730", "-0.73"),
             ("K562", "-0.29"), ("HepG2", "-0.34")):
    add("Fig. 3", "e", f"{k} 5mC ρ", round(f3["meth"][k][0], 2), t, "t9d_coupling", 0.005)

# ---------------- 补充图 ----------------
add("Fig. S1", "a", "预设阈值数", "nine", "nine pre-specified cut-offs", "t1_did_perm4_c5")
add("Fig. S1", "a", "范围外阈值", ", ".join(f"{v:.2f}" for v in s1["outside"]),
    "Cut-offs 0.70 and 0.75 lie outside", "t1_did_perm4_c5")
add("Fig. S1", "b", "C(41,4)", f'{s1["n_comb"]:,}', "C(41,4) = 101,270", "t1_did_perm4_c5")
add("Fig. S1", "b", "p 范围", f'{s1["perm_p_prespecified"][0]:.1e}-{s1["perm_p_prespecified"][1]:.1e}',
    "9.9×10⁻⁶ to 5.9×10⁻⁵", "t1_did_perm4_c5")
add("Fig. S4", "a", "样本数", s2["n_samples"], "13 samples", "t10_sim_band")
add("Fig. S4", "a", "H9 高出带", round(s2["h9_excess"], 3), "0.024 above its band", "t10_sim_band", 0.0005)
add("Fig. S2", "a", "Hap1 双零比例", round(s3["hap1_frac_both0"] * 100, 1), "90.3%", "t15_reliability_41", 0.05)
add("Fig. S2", "b", "GM12878 hap2 覆盖", round(s3["gm12878_nz_hap2"] * 100, 1), "28.0%", "t15_reliability_41", 0.05)
add("Fig. S2", "b", "GM12878 ρ(池化,hap2)", round(s3["gm12878_rho_h2"], 2), "ρ = 0.12", "t15_reliability_41", 0.005)
add("Fig. S2", "b", "39 个样本自洽区间",
    f'{s3["selfcons_range"][0]:.2f}-{s3["selfcons_range"][1]:.2f}', "0.46-0.98", "t15_selfconsistency", 0.005)
add("Fig. S2", "b", "自洽中位", round(s3["selfcons_median"], 2), "median 0.86", "t15_selfconsistency", 0.005)
add("Fig. S5", "a", "n", s4["n"], "946 L1PA2 copies", "t9i_dA_by_copy")
add("Fig. S5", "a", "ρ", round(s4["rho"], 2), "-0.28 (95% CI -0.33 to -0.22", "t9i_dA_by_copy", 0.005)

add("Fig. S5", "b", "五分位序列", " / ".join(f"{v:.2f}" for v in s4["quintiles"]),
    "-0.02, -0.09, -0.21, -0.22 and -0.50", "t9i_dA_by_copy")
add("Fig. S3", "a", "格子数", s5["n_cells"], "31 cells", "t10_readability_map")
add("Fig. S3", "a", "31 格 Spearman", round(s5["rho31"], 2), "Spearman +0.87", "t10_readability_map", 0.005)
add("Fig. S3", "a", "22 格 Spearman", round(s5["rho_sub"], 2), "+0.77", "t10_readability_map", 0.005)
add("Fig. S3", "a", "家族分布", str(s5["by_family"]), "31", "t10_readability_map")

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "t16_figure_number_audit.tsv"), sep="\t", index=False)
bad = df[df.status != "OK"]
print(df[["figure", "panel", "quantity", "value_on_figure", "value_in_manuscript", "status"]]
      .to_string(index=False, max_colwidth=48)[:7000])
print(f"\n共 {len(df)} 项，MISMATCH {len(bad)}")
if len(bad):
    print(bad[["figure", "panel", "quantity", "value_on_figure", "value_in_manuscript", "note"]]
          .to_string(index=False))
