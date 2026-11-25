#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
375 — Table S4 的覆盖感知版本（r26）。

在 t15_table_s3.tsv（41 行）后追加四列：

  n_cov              两条 hap 轨道在该窗都有覆盖的 L1PA2 拷贝数
  frac_both0_covered 双零拷贝中两条 hap 都有覆盖的比例（有覆盖 = 真零，不是无定相）
  r_s_cov            只保留 n_cov 这些拷贝时的 Spearman r_s
  r_full_cov         r_s_cov 的 Spearman–Brown 校正值

只有 26 个样本有 hap 覆盖轨道（25 个可评估样本 + GM12878），其余 15 行为空。
输出：03_results/t15/t15_table_s4_cov.tsv
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")

pa = pd.read_csv(os.path.join(T15, "t15_pa5_hap_L1PA2_41.tsv.gz"), sep="\t").set_index("unit_id")
cv = pd.read_csv(os.path.join(T15, "t15_hap_cov_L1PA2_41.tsv.gz"), sep="\t").set_index("unit_id")
s4 = pd.read_csv(os.path.join(T15, "t15_table_s3.tsv"), sep="\t")

rows = {}
tot_both0 = tot_cov = 0
for s in [c[:-5] for c in pa.columns if c.endswith(".hap1")]:
    a, b = pa[s + ".hap1"], pa[s + ".hap2"]
    ca, cb = cv[s + ".hap1"], cv[s + ".hap2"]
    if ca.isna().all():
        continue
    both0 = (a == 0) & (b == 0)
    covered = (ca > 0) & (cb > 0)
    n_both0 = int(both0.sum())
    n_cov_both0 = int((both0 & covered).sum())
    keep = covered
    rs = spearmanr(a[keep], b[keep]).statistic
    rf = max(0.0, 2 * rs / (1 + rs)) if rs > -1 else 0.0
    rows[s] = dict(n_cov=int(keep.sum()),
                   frac_both0_covered=round(n_cov_both0 / n_both0, 4) if n_both0 else np.nan,
                   r_s_cov=round(float(rs), 4),
                   r_full_cov=round(float(rf), 4),
                   _n_both0=n_both0, _n_cov_both0=n_cov_both0)
    tot_both0 += n_both0
    tot_cov += n_cov_both0

add = pd.DataFrame(rows).T.rename_axis("sample").reset_index()
out = s4.merge(add.drop(columns=["_n_both0", "_n_cov_both0"]), on="sample", how="left")
cols = list(s4.columns) + ["n_cov", "frac_both0_covered", "r_s_cov", "r_full_cov"]
out = out[cols]
out.to_csv(os.path.join(T15, "t15_table_s4_cov.tsv"), sep="\t", index=False)

# 汇总行，供 345 的文字审计直接核对正文里的 56% / 28-92% / 0.00-0.21
assess0 = out[out.n_cov.notna()]
summary = pd.DataFrame([dict(
    n_with_track=len(assess0),
    n_assessable=int(assess0.phased.sum()),
    n_both0=tot_both0,
    n_both0_covered=tot_cov,
    pct_both0_covered=round(100 * tot_cov / tot_both0, 2),
    frac_min=float(assess0[assess0.phased].frac_both0_covered.min()),
    frac_max=float(assess0[assess0.phased].frac_both0_covered.max()),
    r_s_min=float(assess0[assess0.phased].r_s_cov.min()),
    r_s_max=float(assess0[assess0.phased].r_s_cov.max()),
    n_le0=int((assess0[assess0.phased].r_s_cov <= 0).sum()),
    n_ge_half=int((assess0[assess0.phased].r_s_cov >= 0.5).sum()),
    r_full_max=float(assess0[assess0.phased].r_full_cov.max()),
)])
summary.to_csv(os.path.join(T15, "t15_table_s4_cov_summary.tsv"), sep="\t", index=False)
print(summary.T.to_string())

assess = out[(out.n_cov.notna()) & (out.phased == True)]  # noqa: E712  (读出来是 bool)
fin = out[out.n_cov.notna()]
print("samples with coverage track:", len(fin), "| assessable:", len(assess))
print("double-zero copies, all with coverage:", tot_both0,
      "| of which covered on both haplotypes:", tot_cov,
      "(%.1f%%)" % (100 * tot_cov / tot_both0))
print("r_s_cov (assessable): %.4f to %.4f" % (assess.r_s_cov.min(), assess.r_s_cov.max()))
print("r_full_cov max: %.4f" % assess.r_full_cov.max())
print("at or below zero:", int((assess.r_s_cov <= 0).sum()),
      "| at or above 0.5:", int((assess.r_s_cov >= 0.5).sum()))
print("frac_both0_covered range: %.2f to %.2f" % (
    assess.frac_both0_covered.min(), assess.frac_both0_covered.max()))
print("readable labels among assessable:", assess.readable.value_counts().to_dict())
print("wrote", os.path.join(T15, "t15_table_s4_cov.tsv"))
