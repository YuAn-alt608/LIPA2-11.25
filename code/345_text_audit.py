#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
345 — 对稿件（默认 r17，可用 argv 指向 r18）做「文字/表格 vs 结果表」的逐项核对（只读，不改稿）。

r17 的 Table 1 / Table 2 是真正的 docx 表格，所以这里既读段落也读单元格。
输出：03_results/t16/t16_text_audit_r17.tsv
"""
import os
import re
import numpy as np
import pandas as pd
import docx

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")
OUT = os.path.join(PAPER, "docs", "audit")
import sys
DOCX = (sys.argv[1] if len(sys.argv) > 1 else
        os.path.join(PAPER, "manuscript",
                     "Reliability-aware_Fiber-seq_L1PA2_manuscript.docx"))
TAG = os.path.splitext(os.path.basename(DOCX))[0].rsplit("_", 1)[-1]
os.makedirs(OUT, exist_ok=True)

d = docx.Document(DOCX)
paras = [p.text.strip() for p in d.paragraphs]
flat = re.sub(r"\s+", " ", " ".join(paras))
T1 = d.tables[0]
T2 = d.tables[1]
ALLTXT = flat + " " + " ".join(c.text.strip() for t in d.tables for r in t.rows for c in r.cells)

rows = []


def chk(item, claim, truth, ok, note=""):
    rows.append(dict(item=item, claim_in_manuscript=str(claim), value_in_source=str(truth),
                     match="OK" if ok else "MISMATCH", note=note))


rel = pd.read_csv(os.path.join(T15, "t15_reliability_41.tsv"), sep="\t",
                  keep_default_na=False, na_values=[""])
rho = pd.read_csv(os.path.join(T15, "t15_rho_vs_h1_41.tsv"), sep="\t").set_index("sample")
rho_cov = pd.read_csv(os.path.join(T15, "t15_rho_vs_h1_cov_41.tsv"), sep="\t").set_index("sample")
rel_cov = pd.read_csv(os.path.join(T15, "t15_reliability_cov_41.tsv"), sep="\t").set_index("sample")
gc = pd.read_csv(os.path.join(RES, "t13", "t13_group_contrast_sample_first.tsv"), sep="\t")
ps = pd.read_csv(os.path.join(RES, "t13", "t13_per_sample_did.tsv"), sep="\t")
acc = pd.read_csv(os.path.join(T15, "t15_pa5_vs_rfull.tsv"), sep="\t").set_index("subset")
amb = pd.read_csv(os.path.join(RES, "t9m", "t9m_partial_rho.tsv"), sep="\t")
coup = pd.read_csv(os.path.join(RES, "t9d", "t9d_coupling.tsv"), sep="\t")

# ---------- 1 计数 ----------
n_all, n_ph = len(rel), int(rel.phased.sum())
n_rd = int(rel.readable.isin(["yes", "borderline"]).sum())
n_no = int((rel.readable == "no").sum())
n_na = int((rel.readable == "n/a").sum())
n_le0 = int(rel.r_s_le0.sum())
chk("总样本 41", "41", n_all, n_all == 41 and "41 samples" in flat)
chk("通过定相 39", "39", n_ph, n_ph == 39 and "39 " in flat)
chk("可读 10", "10", n_rd, n_rd == 10)
chk("不可读 29", "29", n_no, n_no == 29)
chk("未评估 2", "two could not be assessed", n_na, n_na == 2)
chk("r_s_nz≤0 共 14", "14", n_le0, n_le0 == 14)

# ---------- 2 Highlights / Abstract 口径 ----------
chk("摘要 10/29 口径", "10 supported copy-level reading, 29 had measurable signal",
    "10 / 29", "10 supported copy-level reading, 29 had measurable signal" in flat)
chk("摘要 10/29/2", "10 supported ... 29 ... two could not be assessed",
    f"{n_rd}/{n_no}/{n_na}", "and two could not be assessed" in flat)

# ---------- 3 阈值句 ----------
s7 = pd.read_csv(os.path.join(T15, "t15_table_s7_threshold_cov.tsv"), sep="\t")
chk("阈值句分母", "39 samples that pass the validity rule",
    n_ph, "39 samples that pass the validity rule" in flat and n_ph == 39)
chk("阈值句稳定段", "0.55-0.667, 10 readable",
    f"{list(s7[(s7.threshold >= 0.55) & (s7.threshold <= 0.667)].readable)}",
    "from 0.55 to 0.667 with 10 samples readable" in flat
    and set(s7[(s7.threshold >= 0.55) & (s7.threshold <= 0.667)].readable) == {10})
chk("阈值句指向 Table S3", "Table S3", "-", "Table S3" in flat)

# ---------- 4 Table 1 ----------
body1 = [[c.text.strip() for c in r.cells] for r in T1.rows]
hdr1 = body1[0]
bad1 = []
for r in body1[1:]:
    s = r[0]
    if s not in rho.index:
        bad1.append((s, "not in t15")); continue
    q = rho_cov.loc[s]
    for lab, cell, val, tol in (("r_full_cov", r[2], q.r_sub, 6e-4),
                                ("rho_obs", r[3], q.rho_obs, 6e-4)):
        if abs(float(cell) - float(val)) > tol:
            bad1.append((s, lab, cell, val))
    mm = re.match(r"([\d\.]+) \(([\d\.]+)–([\d\.]+)\)", r[4])
    if mm:
        for k, val in zip(mm.groups(), (q.rho_corr, q.lo, q.hi)):
            if abs(float(k) - float(val)) > 6e-4:
                bad1.append((s, "rho_corr", k, val))
    else:
        bad1.append((s, "rho_corr 格式", r[4]))
    if r[5] != rel_cov.loc[s, "readable"]:
        bad1.append((s, "readable", r[5], rel_cov.loc[s, "readable"]))
chk("Table 1 行数", len(body1) - 1, 13, len(body1) - 1 == 13)
chk("Table 1 数值（13 行 × 5 列）", "见表", "t15_reliability_cov_41 / t15_rho_vs_h1_cov_41",
    not bad1, str(bad1[:5]) if bad1 else "|Δ|≤6e-4")

# ---------- 4b S4 / S5 必须与 Table 1 同一估计器 ----------
_s4c = pd.read_csv(os.path.join(T15, "t15_table_s4_cov.tsv"), sep="\t").set_index("sample")
_s5c = pd.read_csv(os.path.join(T15, "t15_table_s8_pairwise_cov.tsv"), sep="\t")
_d4 = [(s, round(float(_s4c.loc[s, "r_full_cov"]), 4), round(float(rho_cov.loc[s, "r_sub"]), 4))
       for s in rho_cov.index
       if s in _s4c.index and np.isfinite(rho_cov.loc[s, "r_sub"])]
_bad4 = [x for x in _d4 if abs(x[1] - x[2]) > 6e-4]
chk("Table S4 的 r_full,cov 与 Table 1 同值", f"{len(_d4)} 行", f"{len(_bad4)} 处不符", not _bad4,
    str(_bad4[:4]) if _bad4 else "|Δ|≤6e-4")
_bad5 = []
for _, r in _s5c.iterrows():
    if "H1" not in (r.sample_a, r.sample_b):
        continue
    o = r.sample_b if r.sample_a == "H1" else r.sample_a
    if o in rho_cov.index and np.isfinite(r.rho_corr) and \
            abs(float(r.rho_corr) - float(rho_cov.loc[o, "rho_corr"])) > 6e-4:
        _bad5.append((o, round(float(r.rho_corr), 4), round(float(rho_cov.loc[o, "rho_corr"]), 4)))
chk("Table S5 的 H1 列与 Table 1 同值", "10 对", f"{len(_bad5)} 处不符", not _bad5,
    str(_bad5[:4]) if _bad5 else "|Δ|≤6e-4")
chk("Table S5 有 55 对且 ρ_corr 36 对", "55 / 36", f"{len(_s5c)} / {int(_s5c.rho_corr.notna().sum())}",
    len(_s5c) == 55 and int(_s5c.rho_corr.notna().sum()) == 36)

# ---------- 5 Table 2 ----------
body2 = [[c.text.strip() for c in r.cells] for r in T2.rows]
bad2 = []
for r in body2[1:]:
    s = r[0]
    for col, scope in ((1, "ALL"), (2, "INTERGENIC")):
        q = amb[(amb.scope == scope) & (amb.pred == s)]
        if not len(q):
            bad2.append((s, scope, "missing")); continue
        q = q.iloc[0]
        mm = re.match(r"([+\-−]?[\d\.]+) \(([+\-−]?[\d\.]+)–([+\-−]?[\d\.]+)\)", r[col])
        if not mm:
            bad2.append((s, scope, "格式", r[col])); continue
        vals = [float(x.replace("−", "-")) for x in mm.groups()]
        for got, want in zip(vals, (q.rho_partial, q.lo, q.hi)):
            if abs(got - want) > 6e-4:
                bad2.append((s, scope, got, want))
chk("Table 2 行数", len(body2) - 1, 8, len(body2) - 1 == 8)
chk("Table 2 数值（8 行 × 2 列 × 3 数）", "见表", "t9m_partial_rho", not bad2,
    str(bad2[:5]) if bad2 else "|Δ|≤6e-4")
chk("Table 2 分母", "n = 617 / 403", (amb[amb.scope == "ALL"].n.iloc[0], amb[amb.scope == "INTERGENIC"].n.iloc[0]),
    "n = 617" in body2[0][1] and "n = 403" in body2[0][2])

# ---------- 6 组间对比 ----------
for fam, key, stat in (("L1PA2", "5p1k", "raw"), ("L1HS", "5p1k", "did"), ("L1PA3", "5p1k", "did"),
                       ("LTR7", "start400", "raw"), ("AluS", "start400", "raw")):
    r = gc[(gc.window == key) & (gc.family == fam)].iloc[0]
    v = r.raw_diff if stat == "raw" else r.did
    chk(f"{fam} {stat}", f"{v:.1f}", f"{v:.4f}", f"{v:.1f}" in flat,
        f"{r.raw_lo:.1f}–{r.raw_hi:.1f}")
r = gc[(gc.window == "5p1k") & (gc.family == "L1PA2")].iloc[0]
chk("L1PA2 组均值 16.0 / 3.4", "16.0 ... 3.4", f"{r.group_pluri:.2f} / {r.group_other:.2f}",
    "16.0" in flat and "3.4" in flat)
_didr = pd.read_csv(os.path.join(T15, "t15_did_r28.tsv"), sep="\t")
_didcov = _didr.cov_median[_didr.grp == 'pluripotent'].mean() - _didr.cov_median[_didr.grp == 'other'].mean()
_didsub = _didr.raw_median_subfam[_didr.grp == 'pluripotent'].mean() - _didr.raw_median_subfam[_didr.grp == 'other'].mean()
_didfix = _didr.raw_median[_didr.grp == 'pluripotent'].mean() - _didr.raw_median[_didr.grp == 'other'].mean()
chk("L1PA2 DiD 12.3（合并参照）", "+12.3 percentage points", f"{r.did:.2f}",
    "+12.3 percentage points" in flat and "no copy-level interval" in flat)
chk("DiD 合并参照 12.27", "12.266 fixed", f"{_didfix:.3f}",
    abs(_didfix - 12.266) < 0.005 and abs(_didsub - 12.144) < 0.005)
chk("分母口径差 <0.01", "12.265 covered", f"{_didcov:.3f}",
    abs(_didcov - _didfix) < 0.01 and "12.265 covered" in flat)

# ---------- 7 逐样本与留一 ----------
diffs = {r["sample"]: r["diff"] for _, r in ps.iterrows()}
chk("逐样本范围 8.76-15.92", "8.76 ... 15.92",
    "/".join(f"{diffs[s]:.2f}" for s in ("H1", "WTC-11", "H9", "GM27730")),
    "8.76" in flat and "15.92" in flat
    and abs(diffs["GM27730"] - 8.76) < 5e-3 and abs(diffs["H1"] - 15.92) < 5e-3,
    "中间两值 14.47/13.09 只进 Fig. 1e 源数据，不在正文")
chk("逐样本 14.47 / 13.09 仍在源数据", "Fig. 1e",
    f"{diffs['WTC-11']:.2f}/{diffs['H9']:.2f}",
    abs(diffs["WTC-11"] - 14.47) < 5e-3 and abs(diffs["H9"] - 13.09) < 5e-3)
loo = ps.did_after_loo.dropna()
chk("留一 11.3–13.7", "11.3–13.7", f"{loo.min():.1f}–{loo.max():.1f}",
    f"{loo.min():.1f}–{loo.max():.1f}" in flat)

# ---------- 8 可靠性正文 ----------
h = rel.set_index("sample")
for s, a, b in (("Hap1", 0.790, 0.121), ("Panc1", 0.691, 0.258), ("A549", 0.647, 0.300),
                ("THP-1", 0.622, 0.302)):
    ok = abs(h.loc[s, "r_s_all"] - a) < 2e-3 and abs(h.loc[s, "r_s_nz"] - b) < 2e-3
    chk(f"{s} r_s {a}→{b}", f"{a}/{b}",
        f"{h.loc[s,'r_s_all']:.3f}/{h.loc[s,'r_s_nz']:.3f}", ok)
lcl = rel[rel.group == "other (LCL)"]
_lcl_cov = pd.read_csv(os.path.join(T15, "t15_table_s4_cov.tsv"), sep="\t")
_lcl_cov = _lcl_cov[(_lcl_cov.group == "other (LCL)") & (_lcl_cov.phased == True)]  # noqa: E712
chk("25 个 LCL 的 r_s_cov ≤0.14", "at most 0.14",
    f"{_lcl_cov.r_s_cov.max():.3f}",
    abs(_lcl_cov.r_s_cov.max() - 0.14) <= 0.005
    and "reaches at most 0.14 under coverage" in flat)
# r29：正文改报规则实际检验的量（两 hap 均无覆盖），不再报双零
_s4h = pd.read_csv(os.path.join(T15, "t15_table_s4_cov.tsv"), sep="\t").set_index("sample")
chk("Hap1 无覆盖 88.3%", "88.3%", f"{_s4h.loc['Hap1', 'frac_nocov']*100:.1f}%",
    abs(float(_s4h.loc["Hap1", "frac_nocov"]) * 100 - 88.3) < 0.05)
chk("GM12878 hap2 28.0%", "28.0%", f"{h.loc['GM12878','nz_hap2']*100:.1f}%",
    abs(h.loc["GM12878", "nz_hap2"] * 100 - 28.0) < 0.05)
chk("GM12878 ρ(pooled,hap2)=0.12", "0.12", f"{h.loc['GM12878','rho_m_h2']:.2f}",
    abs(h.loc["GM12878", "rho_m_h2"] - 0.12) < 0.005)
_prc = pd.read_csv(os.path.join(T15, "t15_pa5_vs_rfull_cov.tsv"), sep="\t").set_index("subset")
for k, want in (("all assessed", "+0.70"), ("non-pluripotent", "+0.58"), ("LCL", "+0.22")):
    chk(f"ρ(pa5,r_full,cov) {k}", want, f"{_prc.loc[k,'rho']:+.2f}", f"{want}" in flat)
chk("LCL 的 +0.22 与 p 值成句", "+0.22 (p = 0.30)", "text",
    f"+0.22 (p = {_prc.loc['LCL','p']:.2f})" in flat)

# ---------- 9 Fig 3 相关 ----------
for s, want in (("H1", "+0.441"), ("K562", "+0.067"), ("HG002", "−0.020")):
    q = amb[(amb.scope == "ALL") & (amb.pred == s)].iloc[0]
    chk(f"Table 2 {s} all", want, f"{q.rho_partial:+.3f}", f"{want}" in ALLTXT)
chk("WGBS −0.23", "−0.23", "text", "−0.23" in flat)
# H1 这一句在 r18 里由 −0.83(−0.85 to −0.80) 改为 −0.82(−0.84 to −0.80)，用整句匹配防漏
chk("5mC vs H1 整句", "correlated at ρ = −0.82 (95% CI −0.84 to −0.80)", "text",
    "correlated at ρ = −0.82 (95% CI −0.84 to −0.80)" in flat)
for s, want in (("H9", "−0.82"), ("WTC-11", "−0.80"), ("GM27730", "−0.73"),
                ("K562", "−0.29"), ("HepG2", "−0.34")):
    q = coup[coup.pa5_of == s].iloc[0]
    chk(f"5mC vs {s}", want, f"{q.rho:.2f}", f"{want}" in flat)

# ---------- 10 摘要/讨论的派生数 ----------
chk("17 intergenic positives", "17", 17, "17 intergenic positives" in flat or "17 of 403" in flat)
chk("21.8% zero", "21.8%", 21.8, "21.8%" in flat)
chk("315 detected", "315", 315, "315 detected" in flat)
chk("946 copies", "946", 946, "946" in flat)
chk("58.4% 反号", "58.4%", 58.4, "58.4%" in flat)
chk("quintiles", "−0.02, −0.09, −0.21, −0.22, −0.50", "-",
    "−0.02, −0.09, −0.21, −0.22 and −0.50" in flat)

# ---------- 11 Data integrity ----------
chk("246 窗 7.1e-15", "246 windows, maximum difference 7.1×10⁻¹⁵", "-",
    "246 windows across families" in flat and "7.1×10⁻¹⁵" in flat)
chk("排查日记已移出正文", "no debug narrative", "-",
    "The same error recurred" not in flat and "stride equal to the thread count" not in flat)
chk("行序问题只留一句", "permuted row order", "-",
    "permuted row order across files" in flat)

# ---------- 12 覆盖感知主估计量（r28） ----------
cov = pd.read_csv(os.path.join(T15, "t15_table_s4_cov_summary.tsv"), sep="\t").iloc[0]
s4c = pd.read_csv(os.path.join(T15, "t15_table_s4_cov.tsv"), sep="\t")
for _c in ("r_s_cov", "r_full_cov", "frac_both0_covered", "med_frac_cov"):
    s4c[_c] = pd.to_numeric(s4c[_c], errors="coerce")
assess = s4c[s4c.n_cov.notna() & (s4c.phased == True)]  # noqa: E712
notread = assess[assess.readable == "no"]
readable_ = assess[assess.readable.isin(["yes", "borderline"])]
lcl = assess[assess.group == "other (LCL)"]

chk("Table S4 41 行有 n_cov", "41 rows with n_cov filled", int(cov.n_with_track),
    int(cov.n_with_track) == 41, "381 补的 15 个样本 + 331 的 26 个")
chk("可评估 39", "39 assessable samples", int(cov.n_assessable),
    int(cov.n_assessable) == 39 and "39 assessable samples" in flat)
chk("双零有覆盖 27.4%", "27.4%", f"{cov.pct_both0_covered:.1f}%",
    f"{cov.pct_both0_covered:.1f}% are covered on both haplotypes" in flat)
chk("单样本 2.0–92%，中位 40.5%", "2.0% to 92%, median 40.5%",
    f"{cov.frac_min*100:.1f}-{cov.frac_max*100:.0f}%, med {cov.frac_median*100:.1f}%",
    (f"from {cov.frac_min*100:.1f}% to {cov.frac_max*100:.0f}%" in flat
     and f"median {cov.frac_median*100:.1f}%" in flat))
chk("双零覆盖率与可及性 −0.66", "Spearman −0.66", "-",
    "Spearman −0.66 over the 40 samples" in flat)
chk("不可读 29 个 r_full_cov 0.003–0.51", "0.003–0.51",
    f"{notread.r_full_cov.min():.3f}-{notread.r_full_cov.max():.2f}",
    (f"{notread.r_full_cov.min():.3f}–{notread.r_full_cov.max():.2f}" in flat
     and len(notread) == 29))
chk("可读 10 个 r_full_cov 0.68–0.94", "0.68–0.94",
    f"{readable_.r_full_cov.min():.2f}-{readable_.r_full_cov.max():.2f}",
    (f"{readable_.r_full_cov.min():.2f}–{readable_.r_full_cov.max():.2f}" in flat
     and len(readable_) == 10))
chk("无一 r_s_cov ≤0", "none at or below zero", int(cov.n_le0),
    int(cov.n_le0) == 0 and "at or below zero" in flat)
chk("4 个多能样本双零覆盖 2.5–30%", "2.5–30%", "-",
    "cover 2.5–30% of their double zeros" in flat)
chk("LCL 双零覆盖达 92%", "up to 92%", "-",
    "the lymphoblastoid lines up to 92%" in flat)
# r31：LCL 26 个（25 个通过 + GM12878 不通过），正文与 Fig. S1 题注都用 26
n_lcl = int((s4c.group == "other (LCL)").sum())   # 全表 26 个，含未评估的 GM12878
chk("LCL 共 26 个（25 通过）", "26 of the 41",
    f"{n_lcl} LCL / {int(lcl.phased.sum())} pass",
    n_lcl == 26 and int(lcl.phased.sum()) == 25
    and flat.count("26 of the 41 samples are lymphoblastoid lines") == 2
    and "25 of the 41 samples are lymphoblastoid lines" not in flat)
# r31：全文不再把差别归给测序深度（没有降采样分析）
chk("无 sequencing depth 归因", "0 处", flat.count("sequencing depth"), "sequencing depth" not in flat)
chk("Discussion 改说 cut-off 的选择", "choice of cut-off", "-",
    "mixes accessibility with the choice of cut-off" in flat)
chk("Introduction 问题一改问 cut-off", "change of cut-off", "-",
    "Which readout survives a change of cut-off?" in flat)
# r31：删掉 therefore；双零方向改为「抬高 -> 去掉后落下」
chk("噪声模拟段无 therefore", "0 处", flat.count("therefore"),
    "is therefore a gap in dynamic range" not in flat)
chk("双零方向：抬高三个样本", "carry ... over the readable line", "-",
    "they carry A549, Panc1 and THP-1 over the readable line" in flat
    and "pushed three samples below the readable line" not in flat)
chk("可及性相关 +0.70 / +0.58 / +0.22", "+0.70 / +0.58 / +0.22", "-",
    "Spearman +0.70" in flat and "+0.58" in flat and "+0.22" in flat)
chk("同供体对照 0.11 / 0.25", "0.11 / 0.25", "-",
    "\u03c1_obs = 0.11, or 0.25" in flat)
# 多能样本的两个区间：0.93–0.98 只属于「相对 H1」，六对样本是 0.92–0.98
_s5p = pd.read_csv(os.path.join(T15, "t15_table_s8_pairwise_cov.tsv"), sep="\t")
_pl = ["H1", "H9", "WTC-11", "GM27730"]
_six = _s5p[_s5p.sample_a.isin(_pl) & _s5p.sample_b.isin(_pl)].rho_corr
chk("六对多能样本 0.92–0.98", f"{_six.min():.3f}-{_six.max():.3f}",
    f"{_six.min():.2f}–{_six.max():.2f}",
    f"{_six.min():.2f}–{_six.max():.2f}" in flat and len(_six) == 6)
_bare93 = [m.start() for m in re.finditer("0\.93–0\.98", flat)
           if "H1" not in flat[max(0, m.start() - 95):m.start()]]
chk("0.93–0.98 只出现在「相对 H1」处", "0 处裸用", len(_bare93), not _bare93,
    f"{len(_bare93)} 处未标 H1")
chk("DiD 分层 bootstrap", "9.4 to 14.6", "-", "9.4 to 14.6" in flat)
chk("DiD 精确置换 p", "p = 0.003", "-", "p = 0.003" in flat)
chk("DiD 逐样本范围 8.8–15.9", "8.8 to 15.9", "-",
    "8.8 to 15.9" in flat and "−0.1 to 2.1" in flat)
chk("摘要删掉拷贝层 CI", "no copy-level CI", "-",
    "95% CI 11.7" not in flat and "11.7–12.8" not in flat)
chk("Methods 分母约定", "coverage-aware value", "-",
    "coverage-aware value that divides by the covered bases" in flat)
chk("Methods 800 bp 门槛", "at least 800 bp", "-",
    "at least 800 bp of the" in flat)
chk("旧双零口径只作对照", "alternative estimate", "-",
    "alternative estimate" in flat and "double-zero form reports" in flat)
_pcv = pd.read_csv(os.path.join(T15, "t15_hap_cov_L1PA2_41.tsv.gz"), sep="\t", index_col=0)
_v = _pcv.values.astype(float)
_part = float((( _v > 0) & (_v < 800)).sum()) / _v.size
chk("部分覆盖稀少 0.12%", "0.12% of the 84,624",
    f"{_part*100:.2f}% ({int(((_v>0)&(_v<800)).sum())}/{_v.size})",
    f"{_part*100:.2f}% of the 84,624 copy-haplotypes" in flat)
_r80 = s4c.set_index("sample").r_s_cov
_r1 = pd.to_numeric(s4c.set_index("sample").r_s_cov_1bp, errors="coerce")
_d = float((_r80 - _r1).abs().max())
chk("1 bp 与 800 bp 差 ≤0.01", "at most 0.01 in r_s", f"max |diff| {_d:.4f}",
    _d <= 0.0105 and "at most 0.01 in r_s" in flat)
chk("GM12878 hap1 0.97", "0.97", f"{s4c.set_index('sample').loc['GM12878','rho_m_h1']:.2f}",
    "while haplotype 1 does (0.97)" in flat)
chk("Hap1 88.3% 无覆盖", "88.3%",
    f"{s4c.set_index('sample').loc['Hap1','frac_nocov']*100:.1f}%",
    "88.3% of its copies have neither haplotype covered" in flat)
chk("旧 90.3% 双零已全部退出正文", "0 处", flat.count("90.3"), "90.3" not in flat)
# r31 修订：两条单倍型轨道各只在 8.1% 的拷贝上非零（不是 10%）
_ha = s4c.set_index("sample").loc["Hap1"]
chk("Hap1 每条轨道 8% 非零", "8% of copies",
    f"{_ha.nz_hap1*100:.1f}% / {_ha.nz_hap2*100:.1f}%",
    "non-zero in 8% of copies" in flat and abs(_ha.nz_hap1 * 100 - 8.1) < 0.2)
# r31 修订：四个不可读的非 LCL 样本是 1.77–3.22，下界 1.8
_o4 = s4c[s4c["sample"].isin(["A549", "Panc1", "THP-1", "HG002"])].med_pa5_cov
chk("四个其他样本 1.8–3.2%", "1.8 and 3.2%",
    f"{_o4.min():.2f}–{_o4.max():.2f}",
    "the four others between 1.8 and 3.2%" in flat
    and abs(_o4.min() - 1.8) < 0.05 and abs(_o4.max() - 3.2) < 0.05)


# ---------- 13 图注与图号 ----------
for k in ("Fig. 1", "Fig. 2", "Fig. 3", "Fig. S1", "Fig. S2", "Fig. S3", "Fig. S4", "Fig. S5"):
    chk(f"图注 {k}", "存在", "-",
        any(re.match(r"^" + re.escape(k) + r"(\s|$)", p) for p in paras))
chk("Fig. 2 覆盖 41", "all 41 samples", "-", "for all 41 samples" in flat)
# 图注都要有句末句号（原来四条补充图注漏了）
_caps = [p for p in paras if re.match(r"^Fig\.\s?S?\d+\s", p)]
_bad = [p[:28] for p in _caps if not p.rstrip().endswith((".", "。"))]
chk("全部图注以句号结尾", "0 条缺句号", len(_bad), not _bad, str(_bad))
# Hap1 的 0.79→0.121 不在 Fig. 2a 上（它只画成叉），引图处要说明
chk("Hap1 引图已说明只画叉", "Hap1 appears there only as a cross", "-",
    "Hap1 appears there only as a cross, because it is not assessable" in flat)
chk("Fig. 2b 39", "39 samples that pass the validity rule", "-",
    "r_full,cov with 95% intervals for the 39 samples" in flat)
chk("Fig. 3a 与正文一致（+0.70）", "+0.70", "-", "pass the validity rule (Spearman +0.70)" in flat)
chk("Fig. S3 是 31 格", "31 cells", "-", "31 cells with a valid r_s_nz" in flat)
chk("31 格正文改指 S3", "Fig. S3", "-", "Fig. S3" in flat and "that is what Fig. 3a shows" not in flat)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, f"t16_text_audit_{TAG}.tsv"), sep="\t", index=False)
bad = df[df.match != "OK"]
print(df.to_string(index=False)[:6000])
print(f"\n共 {len(df)} 项，MISMATCH {len(bad)}")
if len(bad):
    print(bad.to_string(index=False))
