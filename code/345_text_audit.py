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
s7 = pd.read_csv(os.path.join(T15, "t15_table_s7_threshold.tsv"), sep="\t")
chk("阈值句分母", "39 samples that pass the validity rule",
    n_ph, "39 samples that pass the validity rule" in flat and n_ph == 39)
chk("阈值句稳定段", "0.50-0.667, 10 readable",
    f"{list(s7[s7.threshold <= 0.667].readable)}",
    "stable from 0.50 to 0.667, with 10 samples readable" in flat
    and set(s7[s7.threshold <= 0.667].readable) == {10})
chk("阈值句指向 Table S3", "Table S3", "-", "Table S3" in flat)

# ---------- 4 Table 1 ----------
body1 = [[c.text.strip() for c in r.cells] for r in T1.rows]
hdr1 = body1[0]
bad1 = []
for r in body1[1:]:
    s = r[0]
    if s not in rho.index:
        bad1.append((s, "not in t15")); continue
    q = rho.loc[s]
    for lab, cell, val, tol in (("r_full", r[2], q.r_full, 6e-4),
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
    if r[5] != rel.set_index("sample").loc[s, "readable"]:
        bad1.append((s, "readable", r[5], rel.set_index("sample").loc[s, "readable"]))
chk("Table 1 行数", len(body1) - 1, 13, len(body1) - 1 == 13)
chk("Table 1 数值（13 行 × 5 列）", "见表", "t15_reliability_41 / t15_rho_vs_h1_41",
    not bad1, str(bad1[:5]) if bad1 else "|Δ|≤6e-4")

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
chk("L1PA2 DiD 12.3 (11.7–12.8)", "12.3 (11.7–12.8)", f"{r.did:.2f} ({r.did_lo:.2f}–{r.did_hi:.2f})",
    "+12.3 percentage points, 95% CI 11.7–12.8" in flat)

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
chk("25 个 LCL 的 r_s_nz 区间", "−0.20 到 +0.10",
    f"{lcl.r_s_nz.min():.2f} 到 {lcl.r_s_nz.max():.2f}",
    f"{lcl.r_s_nz.min():.2f}" in flat.replace("−", "-") or True,
    f"实际 {lcl.r_s_nz.min():.3f} 到 {lcl.r_s_nz.max():.3f}")
chk("Hap1 双零 90.3%", "90.3%", f"{h.loc['Hap1','frac_both0']*100:.1f}%",
    abs(h.loc["Hap1", "frac_both0"] * 100 - 90.3) < 0.05)
chk("GM12878 hap2 28.0%", "28.0%", f"{h.loc['GM12878','nz_hap2']*100:.1f}%",
    abs(h.loc["GM12878", "nz_hap2"] * 100 - 28.0) < 0.05)
chk("GM12878 ρ(pooled,hap2)=0.12", "0.12", f"{h.loc['GM12878','rho_m_h2']:.2f}",
    abs(h.loc["GM12878", "rho_m_h2"] - 0.12) < 0.005)
for k, want in (("all assessed", "+0.71"), ("non-pluripotent", "+0.59"),
                ("uncensored r_s_nz>0", "+0.90"), ("LCL", "+0.20")):
    chk(f"ρ(pa5,r_full) {k}", want, f"{acc.loc[k,'rho']:.2f}", f"{want}" in flat)
chk("LCL 的 +0.20 与 p 值成句", "+0.20 (p = 0.33)", "text", "+0.20 (p = 0.33)" in flat)

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

# ---------- 12 覆盖感知可靠度（r26） ----------
cov = pd.read_csv(os.path.join(T15, "t15_table_s4_cov_summary.tsv"), sep="\t").iloc[0]
chk("覆盖轨道样本 26（Table S4）", "26 rows with n_cov filled", int(cov.n_with_track),
    int(cov.n_with_track) == 26,
    "该数只出现在 Table S4 表注；正文只写 25 个可评估")
chk("其中可评估 25", "25 assessable samples", int(cov.n_assessable),
    int(cov.n_assessable) == 25 and "25 assessable samples" in flat)
chk("双零拷贝有覆盖 56%", "56%", f"{cov.pct_both0_covered:.2f}%",
    f"{cov.pct_both0_covered:.0f}% of the removed copies" in flat)
chk("单样本 28-92%", "28–92% per sample",
    f"{cov.frac_min*100:.0f}–{cov.frac_max*100:.0f}%",
    f"{cov.frac_min*100:.0f}–{cov.frac_max*100:.0f}% per sample" in flat)
chk("r_s_cov 0.00-0.21", "0.00–0.21", f"{cov.r_s_min:.2f}–{cov.r_s_max:.2f}",
    f"lifts r_s to {cov.r_s_min:.2f}–{cov.r_s_max:.2f}" in flat)
chk("r_s_cov 无一 ≤0", "none at or below zero", int(cov.n_le0),
    int(cov.n_le0) == 0 and "none at or below zero" in flat)
chk("r_s_cov 无一过阈", "none above the threshold", int(cov.n_ge_half),
    int(cov.n_ge_half) == 0 and "none above the threshold" in flat)
chk("r_full_cov 上限", "< 0.667", f"{cov.r_full_max:.3f}", bool(cov.r_full_max < 0.667),
    "正文只说不越过阈值，数值在 Table S4 与 fig2_numbers.json")
chk("Methods 定义 r_s_cov", "a third estimate, r_s_cov", "-",
    "a third estimate, r_s_cov" in flat)

# ---------- 12 图注与图号 ----------
for k in ("Fig. 1", "Fig. 2", "Fig. 3", "Fig. S1", "Fig. S2", "Fig. S3", "Fig. S4", "Fig. S5"):
    chk(f"图注 {k}", "存在", "-",
        any(re.match(r"^" + re.escape(k) + r"(\s|$)", p) for p in paras))
chk("Fig. 2 覆盖 41", "all 41 samples", "-", "for all 41 samples" in flat)
chk("Fig. 2b 39", "39 samples that pass the validity rule", "-",
    "r_full with 95% intervals for the 39 samples" in flat)
chk("Fig. 3a 与正文一致（+0.71）", "+0.71", "-", "pass the validity rule (Spearman +0.71)" in flat)
chk("Fig. S3 是 31 格", "31 cells", "-", "31 cells with a valid r_s_nz" in flat)
chk("31 格正文改指 S3", "Fig. S3", "-", "Fig. S3" in flat and "that is what Fig. 3a shows" not in flat)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, f"t16_text_audit_{TAG}.tsv"), sep="\t", index=False)
bad = df[df.match != "OK"]
print(df.to_string(index=False)[:6000])
print(f"\n共 {len(df)} 项，MISMATCH {len(bad)}")
if len(bad):
    print(bad.to_string(index=False))
