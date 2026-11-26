#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
398 — Table S3（阈值敏感性）改用 r_full,cov 重算。

旧的 t15_table_s7_threshold.tsv 由 333 用双零口径的 r_full 生成（11 月 25 日）。
图与正文已换到 r_full,cov，因此阈值表也要重出，否则"A549 在 0.5 处可读"这类
计数与正文不一致。

输出 t15/t15_table_s7_threshold_cov.tsv，列与旧表一致。
"""
import os
import numpy as np
import pandas as pd

T15 = "../03_results/t15"
REL = os.path.join(T15, "t15_reliability_cov_41.tsv")
OUT = os.path.join(T15, "t15_table_s7_threshold_cov.tsv")
THR = [0.50, 0.55, 0.60, 0.65, 0.667, 0.70, 0.75, 0.80, 0.90]


def main():
    d = pd.read_csv(REL, sep="\t")
    v = d[d["valid"].astype(str).str.lower().isin(("true", "1"))].copy()
    print(f"可评估样本 {len(v)}")
    rows = []
    for t in THR:
        ok = v["r_full_cov"] >= t
        above = v["r_full_cov_lo"] >= t
        below = v["r_full_cov_hi"] < t
        span = (~above) & (~below)
        rows.append(dict(
            threshold=t,
            n_phased=len(v),
            readable=int(ok.sum()),
            ci_entirely_above=int(above.sum()),
            ci_spans_threshold=int(span.sum()),
            below_point=int((~ok).sum()),
            ci_entirely_below=int(below.sum()),
            readable_samples=", ".join(sorted(v.loc[ok, "sample"])),
            borderline_samples=", ".join(sorted(v.loc[span, "sample"]))))
    s = pd.DataFrame(rows)
    s.to_csv(OUT, sep="\t", index=False)
    print(s.to_string(index=False))
    print("\n→", OUT)


if __name__ == "__main__":
    main()
