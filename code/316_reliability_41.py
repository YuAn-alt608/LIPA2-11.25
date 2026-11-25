#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
316 — 41 个样本的拷贝级可靠性（把 Table 1 / Fig. 2 从 15 扩到 41）。

输入
----
t9b/t9b_pa5_L1PA2_41_fixed.tsv.gz    41 样本合并轨道 pa5（已有，1,032 拷贝）
t13b/t13b_pa5_hap_L1PA2_41.tsv.gz    41 样本 hap1/hap2 pa5（315 产出）

口径（全部沿用 303，一项不改，便于与旧的 15 样本数字对照）
----------------------------------------------------------
r_s_all = ρ(hap1, hap2) 全部拷贝
r_s_nz  = ρ(hap1, hap2) 去掉「两轨同时为 0」的拷贝（主口径）
frac_both0、ρ(merged,hap1)、ρ(merged,hap2) 用于有效性规则（先写死）：
  frac_both0 > 0.5 或 min(ρ(merged,hap)) < 0.2  -> phase_missing
  r_s_nz <= 0                                   -> phase_conflict
r_full = 2 r_s_nz / (1 + r_s_nz)（Spearman–Brown）
可读阈值 r_full >= 0.667；区间跨 0.667 记 borderline
r_full 区间与 ρ_corr 区间：23 条主染色体块 bootstrap，B=1000，两集合各自重采样，
每个重复内重算 r_s_nz 与 r_full（不在重复外拿点估计去校正）。

输出
----
t13b/t13b_reliability_41.tsv     每样本一行（含区间与判读）
t13b/t13b_rho_vs_h1_41.tsv       每样本 ρ_obs / ρ_corr 与区间
t13b/t13b_threshold_41.tsv       阈值 0.5–0.9 的可读计数（全部 41 / 有 hap 分析）
"""
import os
import sys
import numpy as np
import pandas as pd
from importlib import import_module
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m303 = import_module("303_reliability_doublezero")

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
OUT = os.path.join(RES, "t13b")
MERGED = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
HAP = os.path.join(OUT, "t13b_pa5_hap_L1PA2_41.tsv.gz")
UNITS = os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz")
MAN = "../../01_data/meta/fiberseq_samples.tsv"
B = 1000
SEED = 20261124   # 与 303 同种子、同「共同随机数」设计，便于逐位对齐已发表区间
PLURI = ["H1", "H9", "WTC-11", "GM27730"]
THR = [0.5, 0.55, 0.6, 0.65, 0.667, 0.7, 0.8, 0.9]


def sb(r):
    return 2.0 * r / (1.0 + r) if np.isfinite(r) and r > -1 else np.nan


def main():
    os.makedirs(OUT, exist_ok=True)
    hap = pd.read_csv(HAP, sep="\t")
    mrg = pd.read_csv(MERGED, sep="\t")
    units = pd.read_csv(UNITS, sep="\t")
    units = units[units.part == "L1PA_5p"].copy()
    units["grp"] = units.unit_id.str.rsplit("_", n=1).str[-1]
    units = units[units.grp == "L1PA2"].set_index("unit_id").loc[hap.unit_id].reset_index()
    assert list(units.unit_id) == list(hap.unit_id) == list(mrg.unit_id)
    chrom = units.chrom.values

    man = pd.read_csv(MAN, sep="\t").set_index("sample")
    samples = [s for s in man.index]
    assert len(samples) == 41

    # 点估计
    rows = []
    for s in samples:
        a1 = hap[f"{s}.hap1"].values
        a2 = hap[f"{s}.hap2"].values
        m = mrg[s].values
        st = m303.hap_stats(a1, a2, m)
        ok, note = m303.phase_valid(st)
        st.update(sample=s, group=("pluripotent" if s in PLURI else
                                  ("other (LCL)" if str(man.loc[s, "tissue"]).startswith("Lymphoblastoid")
                                   else "other")),
                  valid=bool(ok), note=note)
        st["r_full"] = sb(st["r_s_nz"]) if ok else np.nan
        st["nz_hap1"] = float(np.mean(a1 > 0))
        st["nz_hap2"] = float(np.mean(a2 > 0))
        rows.append(st)
    df = pd.DataFrame(rows)

    # 块 bootstrap
    rng = np.random.default_rng(SEED)
    cats = pd.Categorical(chrom)
    gix = [np.where(cats.codes == j)[0] for j in range(len(cats.categories))]
    main_idx = [i for i, c in enumerate(cats.categories) if str(c) in
                {f"chr{k}" for k in range(1, 23)} | {"chrX"}]
    assert len(main_idx) == 23, len(main_idx)
    # 与 303 一致：同一组块重采样用于所有样本（共同随机数）
    boots = [np.concatenate([gix[j] for j in rng.integers(0, len(gix), len(gix))])
             for _ in range(B)]

    def rs(a1, a2, ix):
        b0 = (a1[ix] == 0) & (a2[ix] == 0)
        j = ix[~b0]
        if len(j) < 30 or np.ptp(a1[j]) == 0 or np.ptp(a2[j]) == 0:
            return np.nan
        return float(spearmanr(a1[j], a2[j]).statistic)

    a1h, a2h = hap["H1.hap1"].values, hap["H1.hap2"].values
    col = {s: (hap[f"{s}.hap1"].values, hap[f"{s}.hap2"].values, mrg[s].values) for s in samples}
    boot = {s: {"rf": [], "rho": []} for s in samples}
    mh1 = col["H1"][2]
    for b, ix in enumerate(boots):
        rf_h1 = sb(rs(a1h, a2h, ix))
        boot["H1"]["rf"].append(rf_h1)
        for s in samples:
            if s == "H1":
                continue
            a1, a2, m = col[s]
            rf = sb(rs(a1, a2, ix))
            boot[s]["rf"].append(rf)
            # ρ_obs：样本合并轨道 vs H1 合并轨道，同一次重复的同一批块
            if np.ptp(m[ix]) > 0 and np.ptp(mh1[ix]) > 0:
                rho = float(spearmanr(m[ix], mh1[ix]).statistic)
            else:
                rho = np.nan
            boot[s]["rho"].append((rho, rf, rf_h1))
        if (b + 1) % 200 == 0:
            print(f"  bootstrap {b+1}/{B}", flush=True)

    out = []
    for _, r in df.iterrows():
        s = r["sample"]
        lo = hi = np.nan
        if np.isfinite(r["r_full"]):
            v = np.array([x for x in boot[s]["rf"] if np.isfinite(x)])
            if len(v) > 50:
                lo, hi = np.nanpercentile(v, [2.5, 97.5])
        r_full = r["r_full"]
        if not np.isfinite(r_full):
            read = "n/a"
        elif r_full >= 0.667:
            read = "borderline" if (np.isfinite(lo) and lo < 0.667) else "yes"
        else:
            read = "no"
        out.append(dict(sample=s, group=r["group"], n=1032, n_nz=r["n_nz"],
                        frac_both0=round(r["frac_both0"], 4),
                        nz_hap1=round(r["nz_hap1"], 4), nz_hap2=round(r["nz_hap2"], 4),
                        med_pa5=round(r["med_pa5_all"], 3),
                        r_s_all=None if not np.isfinite(r["r_s_all"]) else round(r["r_s_all"], 4),
                        r_s_nz=None if not np.isfinite(r["r_s_nz"]) else round(r["r_s_nz"], 4),
                        rho_m_h1=None if not np.isfinite(r["rho_merged_hap1"]) else round(r["rho_merged_hap1"], 4),
                        rho_m_h2=None if not np.isfinite(r["rho_merged_hap2"]) else round(r["rho_merged_hap2"], 4),
                        valid=r["valid"], note=r["note"],
                        r_full=None if not np.isfinite(r_full) else round(r_full, 4),
                        r_full_lo=None if not np.isfinite(lo) else round(lo, 4),
                        r_full_hi=None if not np.isfinite(hi) else round(hi, 4),
                        readable=read))
    rel = pd.DataFrame(out)
    rel.to_csv(os.path.join(OUT, "t13b_reliability_41.tsv"), sep="\t", index=False)

    # ρ vs H1
    h1m = mrg["H1"].values
    rf1 = float(rel.loc[rel["sample"] == "H1", "r_full"].iloc[0])
    rows = []
    for s in samples:
        if s == "H1":
            continue
        m = mrg[s].values
        rho = float(spearmanr(m, h1m).statistic)
        rf = rel.loc[rel["sample"] == s, "r_full"].iloc[0]
        rho_c = rho / np.sqrt(float(rf) * rf1) if np.isfinite(rf) and rf > 0 else np.nan
        v = boot[s]["rho"]
        v = [(a, b, c) for a, b, c in v if np.isfinite(a) and np.isfinite(b) and np.isfinite(c) and b > 0 and c > 0]
        if len(v) > 50:
            cc = np.array([a / np.sqrt(b * c) for a, b, c in v])
            lo, hi = np.nanpercentile(cc, [2.5, 97.5])
        else:
            lo = hi = np.nan
        rows.append(dict(sample=s, rho_obs=round(rho, 4), r_full=rf,
                         rho_corr=None if not np.isfinite(rho_c) else round(rho_c, 4),
                         lo=None if not np.isfinite(lo) else round(lo, 4),
                         hi=None if not np.isfinite(hi) else round(hi, 4),
                         n_boot=len(v), readable=rel.loc[rel["sample"] == s, "readable"].iloc[0]))
    rho_t = pd.DataFrame(rows)
    rho_t.to_csv(os.path.join(OUT, "t13b_rho_vs_h1_41.tsv"), sep="\t", index=False)

    # 阈值敏感性
    th = []
    for t in THR:
        sub = rel[rel["valid"]]
        th.append(dict(threshold=t, n_valid=len(sub),
                       n_pass=int((sub["r_full"] >= t).sum()),
                       n_pass_all41=int((rel["r_full"] >= t).sum())))
    th = pd.DataFrame(th)
    th.to_csv(os.path.join(OUT, "t13b_threshold_41.tsv"), sep="\t", index=False)

    print(rel.to_string(index=False))
    print(th.to_string(index=False))
    print("\n可读样本:", rel.loc[rel.readable == "yes", "sample"].tolist())
    print("borderline:", rel.loc[rel.readable == "borderline", "sample"].tolist())


if __name__ == "__main__":
    main()
