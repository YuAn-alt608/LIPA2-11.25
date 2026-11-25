#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
356 — 生成 Online Resource 1（Tables S1–S8，单文件 xlsx）。

编号按 r21 的首次引用顺序：
S1 样本清单（含更正后的标签）      <- sample_labels.tsv
S2 拷贝集与窗口定义                  <- 手工定义（数字来自 Methods）
S3 各 r_full 阈值下的可读性          <- t15_table_s7_threshold.tsv
S4 41 样本可靠性图谱                <- t15_table_s3.tsv
S5 11 样本两两一致性（55 对）        <- t15_table_s8_pairwise.tsv
S6 偏相关（AMB 记为未表达）          <- t13_partial_rho_allcopies.tsv
S7 跨样本甲基化相关                  <- t9d_coupling.tsv
S8 等位基因层面按区域/家族          <- t10_dA_dM_{regions,families,quintile}.tsv
"""
import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "manuscript",
                   "Online_Resource_1_Tables_S1-S8.xlsx")

HDR = Font(name="Arial", size=10, bold=True)
TTL = Font(name="Arial", size=11, bold=True)
BODY = Font(name="Arial", size=10)
FILL = PatternFill("solid", fgColor="EDEDED")


def sheet(wb, name, title, blocks):
    """blocks: [(subtitle or None, DataFrame or list-of-lists)]"""
    ws = wb.create_sheet(name)
    ncol = max(len(b[1].columns) if isinstance(b[1], pd.DataFrame) else len(b[1][0])
               for b in blocks)
    ws.cell(row=1, column=1, value=title).font = TTL
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=max(ncol, 2))
    r = 3
    for sub, df in blocks:
        if sub:
            ws.cell(row=r, column=1, value=sub).font = HDR
            r += 1
        if isinstance(df, pd.DataFrame):
            cols = list(df.columns)
            for j, c in enumerate(cols, 1):
                cell = ws.cell(row=r, column=j, value=c)
                cell.font, cell.fill = HDR, FILL
                cell.alignment = Alignment(vertical="center", wrap_text=True)
            for _, row in df.iterrows():
                r += 1
                for j, c in enumerate(cols, 1):
                    v = row[c]
                    if pd.isna(v):
                        v = ""
                    elif hasattr(v, "item"):
                        v = v.item()
                    cell = ws.cell(row=r, column=j, value=v)
                    cell.font = BODY
                    if isinstance(v, float):
                        cell.number_format = "0.0000"
            r += 2
        else:
            for j, c in enumerate(df[0], 1):
                cell = ws.cell(row=r, column=j, value=c)
                cell.font, cell.fill = HDR, FILL
            for row in df[1:]:
                r += 1
                for j, v in enumerate(row, 1):
                    ws.cell(row=r, column=j, value=v).font = BODY
            r += 2
    for j in range(1, ncol + 1):
        width = 12
        for row in ws.iter_rows(min_col=j, max_col=j):
            for c in row:
                if isinstance(c.value, str):
                    width = max(width, min(len(c.value) + 2, 55))
        ws.column_dimensions[get_column_letter(j)].width = width
    ws.freeze_panes = "A3"
    return ws


def main():
    wb = Workbook()
    wb.remove(wb.active)

    # ---------------------------------------------------------------- README
    ws = wb.create_sheet("README")
    ws["A1"] = "Online Resource 1 — Tables S1–S8"
    ws["A1"].font = Font(name="Arial", size=12, bold=True)
    rows = [
        ("Sheet", "Content", "Used in"),
        ("S1", "Sample list with corrected labels (41 samples), Cellosaurus checks and per-sample "
               "FIRE coverage", "Methods, Data and sample labels; Results, readout paragraph"),
        ("S2", "Copy sets, windows and readout definitions", "Methods, Copy sets and windows"),
        ("S3", "Readability under each r_full threshold (39 assessable samples)",
         "Results, reliability paragraph"),
        ("S4", "Reliability atlas for all 41 samples (double-zero and pooled-versus-haplotype "
               "diagnostics)", "Table 1 note; Fig. 2; Fig. S2"),
        ("S5", "Pairwise agreement among the 11 samples of the per-copy panel (55 pairs)",
         "Results, copy ranking"),
        ("S6", "Partial correlations with AMB copies counted as non-expressed",
         "Results, RNA anchor"),
        ("S7", "Cross-sample methylation correlations", "Fig. 3e"),
        ("S8", "Allele-level correlations by region, family and |dA| quintile", "Fig. S5"),
        ("", "", ""),
        ("Source", "All values were written by the analysis scripts from the public tracks; no "
                   "number was entered by hand. Column names are the script column names.",
         ""),
        ("Licence", "CC BY 4.0", ""),
    ]
    for i, row in enumerate(rows, start=3):
        for j, v in enumerate(row, 1):
            c = ws.cell(row=i, column=j, value=v)
            c.font = HDR if i == 3 else BODY
            c.alignment = Alignment(vertical="top", wrap_text=True)
            if i == 3:
                c.fill = FILL
    ws.column_dimensions["A"].width = 10
    ws.column_dimensions["B"].width = 90
    ws.column_dimensions["C"].width = 45

    # ---------------------------------------------------------------- S1
    s1 = pd.read_csv(os.path.join(RES, "sample_labels.tsv"), sep="\t")
    dep = pd.read_csv(os.path.join(RES, "t1_did_by_sample.tsv"), sep="\t")[
        ["sample", "n_peaks", "bp_covered", "genome_fire_fraction", "log_depth"]]
    dep = dep.rename(columns={"n_peaks": "FIRE peaks", "bp_covered": "FIRE bp covered",
                              "genome_fire_fraction": "FIRE fraction of genome",
                              "log_depth": "log10 depth (compendium)"})
    s1d = s1.merge(dep, on="sample", how="left")
    sheet(wb, "S1", "Table S1  Sample list with corrected labels and FIRE coverage (41 samples)", [
        ("Labels as used in this study. GM27730 is listed as lymphoblastoid in the compendium "
         "sample sheet and as iPSC by Cellosaurus; the Cellosaurus class is used here.", s1),
        ("Per-sample FIRE coverage, computed from each sample's all.percent.accessible.bw over the "
         "hg38 autosomes. These columns are the sequencing-depth proxy referred to in Results: "
         "GM27730 has 1.7% of its genome covered by FIRE intervals against 5.0-5.3% in the other "
         "three pluripotent samples.", s1d)])

    # ---------------------------------------------------------------- S2
    fam = [["family", "n_copies", "role", "5' 1 kb window", "400 bp window", "whole element"],
           ["L1PA2", 1032, "target", "yes", "yes", "no"],
           ["L1HS", 332, "comparison", "yes", "yes", "no"],
           ["L1PA3", 1489, "comparison", "yes", "yes", "no"],
           ["L1PA4-7", 4048, "reference set", "yes", "yes", "no"],
           ["LTR7", 2166, "positive control", "no", "yes", "no"],
           ["LTR7B", 597, "control", "no", "yes", "no"],
           ["LTR7Y", 190, "control", "no", "yes", "no"],
           ["LTR5_Hs", 558, "control", "no", "yes", "no"],
           ["AluS", 5000, "near-zero control", "no", "yes", "capped at 5 kb"],
           ["AluY", 5000, "control", "no", "yes", "capped at 5 kb"],
           ["SVA_D", 991, "control", "no", "yes", "no"],
           ["SVA_E", 165, "control", "no", "yes", "no"],
           ["SVA_F", 446, "control", "no", "yes", "no"]]
    win = [["window", "definition"],
           ["5' 1 kb", "first 1 kb from the element's biological start (L1 families only)"],
           ["400 bp", "first 400 bp from the biological start (all families; short elements "
                      "comparable)"],
           ["whole element", "entire element, capped at 5 kb for Alu"]]
    meas = [["measure", "definition"],
            ["pa5", "mean percent-accessible over the 5' 1 kb window (0-100; percentage of "
                    "fibres with a FIRE element at a position)"],
            ["c5", "fraction of the same window covered by FIRE intervals (0-1); needs a "
                   "threshold to call a copy accessible"]]
    sheet(wb, "S2", "Table S2  Copy sets, windows and readout definitions", [
        ("Copy sets (UCSC hg38 RepeatMasker annotation)", fam),
        ("Windows; minus-strand copies are anchored at the last base of the element", win),
        ("Readouts", meas)])

    # ---------------------------------------------------------------- S3
    s3 = pd.read_csv(os.path.join(RES, "t15", "t15_table_s7_threshold.tsv"), sep="\t")
    sheet(wb, "S3", "Table S3  Readability under each r_full threshold (39 assessable samples)", [
        ("readable = r_full at or above the threshold; ci_entirely_above and ci_entirely_below "
         "count samples whose 95% interval lies wholly on one side", s3)])

    # ---------------------------------------------------------------- S4
    s4 = pd.read_csv(os.path.join(RES, "t15", "t15_table_s3.tsv"), sep="\t")
    sheet(wb, "S4", "Table S4  Reliability atlas for all 41 samples", [
        ("All 41 samples. r_s_all retains double-zero copies, r_s_nz removes them; phased = "
         "False marks the two samples that fail the validity rule", s4)])

    # ---------------------------------------------------------------- S5
    s5 = pd.read_csv(os.path.join(RES, "t15", "t15_table_s8_pairwise.tsv"), sep="\t")
    sheet(wb, "S5", "Table S5  Pairwise agreement among the 11 samples of the per-copy panel "
                    "(55 pairs)", [
        ("Attenuation-corrected values only where both samples are assessable", s5)])

    # ---------------------------------------------------------------- S6
    s6 = pd.read_csv(os.path.join(RES, "t13", "t13_partial_rho_allcopies.tsv"), sep="\t")
    sheet(wb, "S6", "Table S6  Partial correlations with AMB copies counted as non-expressed", [
        ("scope = ALL or INTERGENIC; pred = sample supplying pa5; n_cov = covariates retained", s6)])

    # ---------------------------------------------------------------- S7
    s7 = pd.read_csv(os.path.join(RES, "t9d", "t9d_coupling.tsv"), sep="\t")
    sheet(wb, "S7", "Table S7  Cross-sample methylation correlations", [
        ("rho = Spearman correlation between H1 Fiber-seq 5mC and pa5 of the sample named in "
         "pa5_of; same_cells marks the same-read comparison", s7)])

    # ---------------------------------------------------------------- S8
    reg = pd.read_csv(os.path.join(RES, "t9i", "t10_dA_dM_regions.tsv"), sep="\t")
    fam2 = pd.read_csv(os.path.join(RES, "t9i", "t10_dA_dM_families.tsv"), sep="\t")
    q = pd.read_csv(os.path.join(RES, "t9i", "t10_dA_dM_quintile.tsv"), sep="\t")
    sheet(wb, "S8", "Table S8  Allele-level correlations in H1 by region, family and |dA| "
                    "quintile", [
        ("By region (rho = Spearman correlation between dA and dM)", reg),
        ("By family", fam2),
        ("By quintile of |dA|", q)])

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    wb.save(OUT)
    print("wrote", OUT, os.path.getsize(OUT) // 1024, "KB")
    for n in wb.sheetnames:
        print("  ", n)


if __name__ == "__main__":
    main()
