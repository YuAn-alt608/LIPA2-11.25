#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
332 — 用修正行序后的 hap 表重算 41 个样本的可靠性（ROUND15）。

输入
----
t9b/t9b_pa5_L1PA2_41_fixed.tsv.gz     41 样本合并轨道 pa5（未受影响）
t15/t15_pa5_hap_L1PA2_41.tsv.gz       331 重抓并逐位校验后的 hap1/hap2

规则（v2，本轮改的只有一条：负 r_s 不再判为「相冲突 / 无效」）
--------------------------------------------------------------
分层是两件事，混在一起会把「测了、不可读」写成「未评估」：

  第一层 定相是否可用（决定该样本进不进分母）
    frac_nocov > 0.5                        -> phase_missing，not assessed
      （frac_nocov = 两条单倍型都没覆盖的拷贝占比；等价于「测得的是空轨道」，
        与双零 frac_both0 不同：双零里还有一条已覆盖、另一条为零的拷贝）
    min(ρ(merged,hap1), ρ(merged,hap2)) < 0.2 -> phase_missing，not assessed

  第二层 可用者里逐拷贝排序是否可靠（决定 readable）
    r_s_nz = ρ(hap1, hap2) 去掉双零拷贝
    r_full = 2 r / (1 + r)（Spearman–Brown），r 截断在 0：
        r_s_nz <= 0  ->  r_full = 0，readable = no，note = `r_s<=0`
    理由：两半轨道的相关在 0 附近抖动是**零结果**，不是「两条轨道互相矛盾」的证据；
    负的信度估计在数学上不可能（信度是方差比），只能说明估不出正值。
    原始负值仍保留在表里（r_s_nz 列），不抹掉。
    r_full >= 0.667 -> readable；区间下界 < 0.667 -> borderline

区间：23 条主染色体块 bootstrap，B=1000，共同随机数（与 303 同 SEED），
每个重复内重算 r_s_nz 再截断取 r_full；ρ_corr 用同一重复的 r_full。

输出（全部在 t15/ 下）
--------------------
t15_reliability_41.tsv     每样本一行
t15_rho_vs_h1_41.tsv       ρ_obs / ρ_corr 与区间
t15_threshold_41.tsv       阈值 0.5–0.9 的可读计数（分母 40 / 41）
t15_table_s3.tsv           Table S3（41 行）
t15_pa5_vs_rfull.tsv       样本层面 ρ(pa5 中位, r_full)，40 点与分成两组
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
OUT = os.path.join(RES, "t15")
MERGED = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
HAP = os.path.join(OUT, "t15_pa5_hap_L1PA2_41.tsv.gz")
HCOV = os.path.join(OUT, "t15_hap_cov_L1PA2_41.tsv.gz")
UNITS = os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz")
MAN = "../../01_data/meta/fiberseq_samples.tsv"
PUB = os.path.join(RES, "t10_reliability_fig.tsv")
B = 1000
SEED = 20261124
PLURI = ["H1", "H9", "WTC-11", "GM27730"]
THR = [0.5, 0.55, 0.6, 0.65, 0.667, 0.7, 0.8, 0.9]
EPS = 1e-12


def sb_clip(r):
    """Spearman–Brown，r 截断在 [0, 1)。r<=0 返回 0（信度不可能为负）。"""
    if not np.isfinite(r):
        return np.nan
    r = min(max(r, 0.0), 0.999)
    return 2.0 * r / (1.0 + r)


def phase_status(st):
    # 判断的是「这条单倍型在这个窗口里根本没有 fibre」，而不是「两条单倍型的 pa5 都是 0」：
    # 双零还包含「有覆盖但确实关闭」的拷贝，后者是可测的，不能当作缺相位剔除。
    if st["frac_uncovered"] > 0.5:
        return False, "phase_missing: both haplotypes uncovered in >50% of copies"
    lo = min(np.nan_to_num(st["rho_merged_hap1"], nan=-1.0),
             np.nan_to_num(st["rho_merged_hap2"], nan=-1.0))
    if lo < 0.2:
        return False, "phase_missing: rho(merged,hap)<0.2"
    if not np.isfinite(st["r_s_nz"]):
        return False, "r_s_nan"
    return True, ""


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
    samples = list(man.index)
    assert len(samples) == 41

    # 覆盖轨道用来判断「有没有测到」，与「测得是 0」分开（有效性规则见 phase_status）
    hcov = pd.read_csv(HCOV, sep="\t").set_index("unit_id").loc[hap.unit_id]

    # ---------- 点估计 ----------
    rows = []
    for s in samples:
        a1 = hap[f"{s}.hap1"].values
        a2 = hap[f"{s}.hap2"].values
        m = mrg[s].values
        st = m303.hap_stats(a1, a2, m)
        st["frac_uncovered"] = float(np.mean((hcov[f"{s}.hap1"].values == 0) &
                                             (hcov[f"{s}.hap2"].values == 0)))
        ok, note = phase_status(st)
        grp = ("pluripotent" if s in PLURI else
               ("other (LCL)" if str(man.loc[s, "tissue"]).startswith("Lymphoblastoid") else "other"))
        r_full = sb_clip(st["r_s_nz"]) if ok else np.nan
        rows.append(dict(sample=s, group=grp, tissue=str(man.loc[s, "tissue"]),
                         n=st["n"], n_nz=st["n_nz"], frac_both0=st["frac_both0"],
                         nz_hap1=float(np.mean(a1 > 0)), nz_hap2=float(np.mean(a2 > 0)),
                         med_pa5=st["med_pa5_all"],
                         r_s_all=st["r_s_all"], r_s_nz=st["r_s_nz"],
                         rho_m_h1=st["rho_merged_hap1"], rho_m_h2=st["rho_merged_hap2"],
                         phased=bool(ok), note=note, r_full=r_full,
                         r_s_used=(min(max(st["r_s_nz"], 0.0), 0.999) if ok else np.nan)))
    df = pd.DataFrame(rows)

    # ---------- 块 bootstrap ----------
    rng = np.random.default_rng(SEED)
    cats = pd.Categorical(chrom)
    gix = [np.where(cats.codes == j)[0] for j in range(len(cats.categories))]
    assert len(gix) == 23, len(gix)
    boots = [np.concatenate([gix[j] for j in rng.integers(0, len(gix), len(gix))])
             for _ in range(B)]

    def rs(a1, a2, ix):
        b0 = (a1[ix] == 0) & (a2[ix] == 0)
        j = ix[~b0]
        if len(j) < 30 or np.ptp(a1[j]) == 0 or np.ptp(a2[j]) == 0:
            return np.nan
        return float(spearmanr(a1[j], a2[j]).statistic)

    col = {s: (hap[f"{s}.hap1"].values, hap[f"{s}.hap2"].values, mrg[s].values)
           for s in samples}
    boot = {s: {"rf": [], "rho": []} for s in samples}
    a1h, a2h = col["H1"][0], col["H1"][1]
    mh1 = col["H1"][2]
    for b, ix in enumerate(boots):
        rf_h1 = sb_clip(rs(a1h, a2h, ix))
        boot["H1"]["rf"].append(rf_h1)
        for s in samples:
            if s == "H1":
                continue
            a1, a2, m = col[s]
            boot[s]["rf"].append(sb_clip(rs(a1, a2, ix)))
            rho = (float(spearmanr(m[ix], mh1[ix]).statistic)
                   if np.ptp(m[ix]) > 0 and np.ptp(mh1[ix]) > 0 else np.nan)
            boot[s]["rho"].append((rho, boot[s]["rf"][-1], rf_h1))
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
        rf = r["r_full"]
        if not r["phased"]:
            read = "n/a"
        elif rf >= 0.667:
            read = "borderline" if (np.isfinite(lo) and lo < 0.667) else "yes"
        else:
            read = "no"
        out.append(dict(sample=s, group=r["group"], n=int(r["n"]), n_nz=int(r["n_nz"]),
                        frac_both0=round(r["frac_both0"], 4),
                        nz_hap1=round(r["nz_hap1"], 4), nz_hap2=round(r["nz_hap2"], 4),
                        med_pa5=round(r["med_pa5"], 3),
                        r_s_all=None if not np.isfinite(r["r_s_all"]) else round(r["r_s_all"], 4),
                        r_s_nz=None if not np.isfinite(r["r_s_nz"]) else round(r["r_s_nz"], 4),
                        rho_m_h1=None if not np.isfinite(r["rho_m_h1"]) else round(r["rho_m_h1"], 4),
                        rho_m_h2=None if not np.isfinite(r["rho_m_h2"]) else round(r["rho_m_h2"], 4),
                        phased=bool(r["phased"]), note=r["note"],
                        r_full=None if not np.isfinite(rf) else round(rf, 4),
                        r_full_lo=None if not np.isfinite(lo) else round(lo, 4),
                        r_full_hi=None if not np.isfinite(hi) else round(hi, 4),
                        readable=read,
                        r_s_le0=bool(r["phased"] and r["r_s_nz"] <= 0)))
    rel = pd.DataFrame(out)
    rel.to_csv(os.path.join(OUT, "t15_reliability_41.tsv"), sep="\t", index=False)

    # ---------- ρ vs H1 ----------
    h1m = mrg["H1"].values
    rf1 = float(rel.loc[rel["sample"] == "H1", "r_full"].iloc[0])
    rows = []
    for s in samples:
        if s == "H1":
            continue
        m = mrg[s].values
        rho = float(spearmanr(m, h1m).statistic)
        rf = rel.loc[rel["sample"] == s, "r_full"].iloc[0]
        ok = np.isfinite(rf) and rf > 0
        rho_c = rho / np.sqrt(float(rf) * rf1) if ok else np.nan
        # 点估计没有定义时区间也没有意义；rf=0 时分母为 0，重复内的比值会发散
        v = [(a, b, c) for a, b, c in boot[s]["rho"]
             if np.isfinite(a) and np.isfinite(b) and np.isfinite(c) and b > 0 and c > 0] if ok else []
        lo, hi = (np.nanpercentile([a / np.sqrt(b * c) for a, b, c in v], [2.5, 97.5])
                  if len(v) > 50 else (np.nan, np.nan))
        rows.append(dict(sample=s, rho_obs=round(rho, 4),
                         r_full=None if not np.isfinite(rf) else rf,
                         rho_corr=None if not np.isfinite(rho_c) else round(rho_c, 4),
                         lo=None if not np.isfinite(lo) else round(lo, 4),
                         hi=None if not np.isfinite(hi) else round(hi, 4),
                         n_boot=len(v), readable=rel.loc[rel["sample"] == s, "readable"].iloc[0]))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "t15_rho_vs_h1_41.tsv"), sep="\t", index=False)

    # ---------- 阈值敏感性 ----------
    th = []
    for t in THR:
        sub = rel[rel["phased"]]
        th.append(dict(threshold=t, n_phased=len(sub),
                       n_pass=int((sub["r_full"] >= t).sum()),
                       n_pass_all41=int((rel["r_full"] >= t).sum())))
    th = pd.DataFrame(th)
    th.to_csv(os.path.join(OUT, "t15_threshold_41.tsv"), sep="\t", index=False)

    # ---------- Table S3 ----------
    s3 = rel.merge(pd.read_csv(os.path.join(OUT, "t15_rho_vs_h1_41.tsv"), sep="\t"),
                   on="sample", how="left", suffixes=("", "_rho"))
    s3["batch"] = np.where(s3["sample"].isin(
        ["H1", "H9", "WTC-11", "GM27730", "Hek293T", "Caco-2", "HepG2", "K562", "MCF-7",
         "Jurkat", "A549", "Panc1", "HG002", "THP-1", "Hap1"]), "round-13 (15)", "round-15 (26)")
    s3 = s3[["sample", "group", "batch", "n", "n_nz", "frac_both0", "nz_hap1", "nz_hap2",
             "med_pa5", "r_s_all", "r_s_nz", "rho_m_h1", "rho_m_h2", "phased", "note",
             "r_full", "r_full_lo", "r_full_hi", "readable", "rho_obs", "rho_corr", "lo", "hi"]]
    s3.to_csv(os.path.join(OUT, "t15_table_s3.tsv"), sep="\t", index=False, na_rep="",
              float_format="%.4f")

    # ---------- 样本层面：可读性 vs 可及性 ----------
    v40 = rel[rel["phased"]].copy()
    lines = []
    for lab, d in ((f"all assessed", v40),
                   ("pluripotent", v40[v40["group"] == "pluripotent"]),
                   ("non-pluripotent", v40[v40["group"] != "pluripotent"]),
                   ("LCL", v40[v40["group"] == "other (LCL)"]),
                   ("uncensored r_s_nz>0", v40[v40["r_s_nz"] > 0]),
                   ("readable", v40[v40["readable"].isin(["yes", "borderline"])])):
        if len(d) >= 4 and d["r_full"].nunique() > 1:
            r = spearmanr(d["med_pa5"], d["r_full"])
            lines.append(dict(subset=lab, n=len(d), rho=round(float(r.statistic), 4),
                              p=float(f"{r.pvalue:.3g}")))
    pr = pd.DataFrame(lines)
    pr.to_csv(os.path.join(OUT, "t15_pa5_vs_rfull.tsv"), sep="\t", index=False)

    # ---------- 与已发表区间核对（旧 15 个应逐位不变） ----------
    if os.path.exists(PUB):
        pub = pd.read_csv(PUB, sep="\t").set_index("sample")
        chk = []
        for _, r in rel.iterrows():
            s = r["sample"]
            if s in pub.index and np.isfinite(r["r_full"]) and np.isfinite(pub.loc[s, "r_full"]):
                chk.append(dict(sample=s,
                                pub_r_full=pub.loc[s, "r_full"], new_r_full=r["r_full"],
                                pub_lo=pub.loc[s, "r_full_lo"], new_lo=r["r_full_lo"],
                                pub_hi=pub.loc[s, "r_full_hi"], new_hi=r["r_full_hi"]))
        ck = pd.DataFrame(chk)
        if len(ck):
            ck["d_lo"] = (ck["pub_lo"] - ck["new_lo"]).abs()
            ck.to_csv(os.path.join(OUT, "t15_vs_published_intervals.tsv"), sep="\t", index=False)
            print(f"与已发表区间核对：{len(ck)} 个样本，max|Δlo| = {ck['d_lo'].max():.4f}")

    # ---------- 打印 ----------
    show = ["sample", "group", "med_pa5", "frac_both0", "rho_m_h1", "rho_m_h2",
            "r_s_nz", "r_full", "r_full_lo", "r_full_hi", "readable"]
    print("\n" + rel[show].to_string(index=False))
    print("\n" + th.to_string(index=False))
    print("\nphased:", int(rel.phased.sum()), "| readable yes:",
          int((rel.readable == "yes").sum()), "| borderline:",
          int((rel.readable == "borderline").sum()), "| no:",
          int((rel.readable == "no").sum()), "| n/a:", int((rel.readable == "n/a").sum()))
    print("r_s_nz<=0 但 phased 的样本:", rel.loc[rel.r_s_le0, "sample"].tolist())
    print("\n样本层面：\n" + pr.to_string(index=False))


if __name__ == "__main__":
    main()
