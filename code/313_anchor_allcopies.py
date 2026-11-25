#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
313 — 全拷贝（AMB 并入）版本的偏相关，带块 bootstrap 区间（ROUND13 P0-5a）。

306 的全拷贝版只有点估计（+0.461 / +0.497）。本轮补区间：
  - 口径：保留 AMB 拷贝，连续 RNA 值**不置 0**，只在二分类终点里算作非 EXP；
  - 偏 Spearman，协变量 7 个（intergenic 6 个），与 306 完全一致；
  - 23 染色体块 bootstrap，B=1000，每个重复内重算秩残差。
同时输出 EXP/AMB/SIL 的操作定义（阈值、聚合、分配规则），供 Methods 使用。

输出：t13/t13_partial_rho_allcopies.tsv
      t13/t13_label_definition.tsv
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module

m306 = import_module("306_anchor_partial")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
OUT = os.path.join(RES, "t13")
T9M = os.path.join(RES, "t9m")
T9B = os.path.join(RES, "t9b")
B = 1000
SEED = 20261124


def main():
    os.makedirs(OUT, exist_ok=True)
    rna = pd.read_csv(os.path.join(T9M, "t9m_rna_by_copy.tsv.gz"), sep="\t")
    fixed = pd.read_csv(os.path.join(T9B, "t9b_pa5_L1PA2_41_fixed.tsv.gz"), sep="\t")
    rna = rna.drop(columns=[c for c in rna.columns if c.startswith("pa5_")])
    rna = rna.merge(fixed, on="unit_id", how="left", validate="one_to_one")
    for s in m306.PRED:
        rna[f"pa5_{s}"] = rna[s].astype(float)
    rna["rna"] = rna[["senseb_r1", "senseb_r2", "senseb_r3"]].median(axis=1)
    rna["chrom"] = rna.chrom.astype(str)

    rows = []
    for scope, sub in (("ALL", rna), ("INTERGENIC", rna[rna.cat == "intergenic"])):
        d = m306.prep(sub.dropna(subset=["rna"]), drop_amb=False)   # AMB 保留
        covs = [c for c in m306.COV if d[c].astype(float).nunique() > 1]
        C = m306.rank_cov(d, covs)
        gix = m306.blocks(d.chrom.values)
        rng = np.random.default_rng(SEED)
        for p in m306.PRED:
            obs = m306.partial_sp(d[f"pa5_{p}"].to_numpy(float), d.rna.to_numpy(float), C)
            vals = []
            for _ in range(B):
                ix = np.concatenate([gix[j] for j in rng.integers(0, len(gix), len(gix))])
                dd = d.iloc[ix]
                Cb = m306.rank_cov(dd, covs)
                vals.append(m306.partial_sp(dd[f"pa5_{p}"].to_numpy(float),
                                            dd.rna.to_numpy(float), Cb))
            v = np.array(vals, float); v = v[np.isfinite(v)]
            rows.append(dict(scope=scope, pred=p, n=len(d), n_cov=len(covs), n_amb=int((d.label_body == "AMB").sum()),
                             rho_partial=obs, lo=float(np.percentile(v, 2.5)),
                             hi=float(np.percentile(v, 97.5)),
                             rho_raw=float(m306.spearmanr(d[f"pa5_{p}"], d.rna).statistic)))
            print(f"{scope:10s} {p:8s} n={len(d):4d} ρ={obs:+.3f} "
                  f"[{np.percentile(v,2.5):+.3f},{np.percentile(v,97.5):+.3f}]", flush=True)
    P = pd.DataFrame(rows)
    P.to_csv(os.path.join(OUT, "t13_partial_rho_allcopies.tsv"), sep="\t", index=False)

    # 标签定义（可复核的计数）
    S = rna[["senseb_r1", "senseb_r2", "senseb_r3"]].to_numpy(float)
    lab = pd.DataFrame([
        dict(item="normalisation", value="每个重复的窗均值除以该重复全基因组均值（CRG 50-mer 可比对区段）"),
        dict(item="window", value="元件全长（el_start–el_end），按链方向取 sense"),
        dict(item="aggregation", value="三个重复的**中位数**（rna）；均值作敏感性"),
        dict(item="EXP rule", value="中位 ≥ 1.0 且至少 2/3 重复 ≥ 0.5"),
        dict(item="SIL rule", value="3/3 重复 < 0.2"),
        dict(item="AMB rule", value="其余（不满足 EXP 也不满足 SIL）"),
        dict(item="counts ALL", value=f"EXP {int((rna.label_body=='EXP').sum())} / "
                                     f"AMB {int((rna.label_body=='AMB').sum())} / "
                                     f"SIL {int((rna.label_body=='SIL').sum())}"),
        dict(item="counts intergenic", value=f"EXP {int(((rna.cat=='intergenic')&(rna.label_body=='EXP')).sum())} / "
                                             f"AMB {int(((rna.cat=='intergenic')&(rna.label_body=='AMB')).sum())} / "
                                             f"SIL {int(((rna.cat=='intergenic')&(rna.label_body=='SIL')).sum())}"),
        dict(item="AMB_as_neg", value="只在二分类终点里把 AMB 记为 0；连续分析的 RNA 值保持原值，不置 0"),
        dict(item="n used (all copies)", value=f"{len(m306.prep(rna.dropna(subset=['rna']), drop_amb=False))}"),
        dict(item="n used (intergenic)", value=f"{len(m306.prep(rna[rna.cat=='intergenic'].dropna(subset=['rna']), drop_amb=False))}"),
    ])
    lab.to_csv(os.path.join(OUT, "t13_label_definition.tsv"), sep="\t", index=False)
    print(lab.to_string(index=False))


if __name__ == "__main__":
    main()
