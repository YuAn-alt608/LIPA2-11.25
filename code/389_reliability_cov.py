#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
389 — r28 主估计量：覆盖感知 + 覆盖碱基分母的 split-half 可靠性。

与 316 的差别（只有这三处，其余口径不变）
------------------------------------------
1. 分母：pa5 原为 Σv / 1,000（未覆盖碱基记 0）。现改为 pa5_obs = Σv / covered_bp，
   即窗口内真正有 fibre 覆盖的碱基作分母，covered_bp 来自 387 重取的池化覆盖与
   t15 的 hap 覆盖。
2. 保留规则：原为「两 hap 在窗内各覆盖 ≥1 bp」。现为「两 hap 各自覆盖 ≥800 bp」
   （窗口的 80%），避免单条 fibre 决定一个拷贝的值。
3. 点估计与 bootstrap 全部用 r_s_cov；r_s_nz 仍写出，作为敏感性列。

区间：23 条主染色体块 bootstrap，B=1000，共同随机数（与 316 同种子）。
每个重复内重算 r_s_cov 与 r_full_cov，不在重复外拿点估计去校正。

输出
----
t15/t15_reliability_cov_41.tsv   每样本一行（点估计、区间、判读）
t15/t15_rho_vs_h1_cov_41.tsv     ρ_obs / ρ_corr 与区间
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")
OUT = T15
HAP = os.path.join(T15, "t15_pa5_hap_L1PA2_41.tsv.gz")
HCOV = os.path.join(T15, "t15_hap_cov_L1PA2_41.tsv.gz")
PCOV = os.path.join(T15, "t15_pool_cov_L1PA2_41.tsv.gz")
POOL = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
UNITS = os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz")
MAN = "../../01_data/meta/fiberseq_samples.tsv"
B = 1000
SEED = 20261124
MINBP = 800          # 两 hap 各需覆盖的碱基数
PLURI = ["H1", "H9", "WTC-11", "GM27730"]


def sb(r):
    return 2.0 * r / (1.0 + r) if np.isfinite(r) and r > -1 else np.nan


def main():
    hap = pd.read_csv(HAP, sep="\t").set_index("unit_id")
    hcov = pd.read_csv(HCOV, sep="\t").set_index("unit_id")
    pcov = pd.read_csv(PCOV, sep="\t").set_index("unit_id")
    pool = pd.read_csv(POOL, sep="\t").set_index("unit_id")
    assert list(hap.index) == list(hcov.index) == list(pcov.index) == list(pool.index)

    units = pd.read_csv(UNITS, sep="\t")
    units = units[(units.part == "L1PA_5p") & (units.unit_id.str.endswith("_L1PA2"))]
    units = units.drop_duplicates("unit_id").set_index("unit_id").loc[hap.index]
    chrom = units.chrom.values
    man = pd.read_csv(MAN, sep="\t")
    samples = [s for s in pool.columns if s != "unit_id"]

    # ---- 逐拷贝观测值（分母修正） ----
    obs_h = {}
    keep = {}
    for s in samples:
        for h in ("hap1", "hap2"):
            p = hap[f"{s}.{h}"].values.astype(float)
            c = hcov[f"{s}.{h}"].values.astype(float)
            obs_h[(s, h)] = np.where(c > 0, p * 1000.0 / np.maximum(c, 1), np.nan)
        keep[s] = (hcov[f"{s}.hap1"].values >= MINBP) & (hcov[f"{s}.hap2"].values >= MINBP)
    obs_pool = {}
    for s in samples:
        p = pool[s].values.astype(float)
        c = pcov[s].values.astype(float)
        obs_pool[s] = np.where(c > 0, p * 1000.0 / np.maximum(c, 1), np.nan)

    def rs_cov(s, ix):
        k = keep[s][ix] & np.isfinite(obs_h[(s, "hap1")][ix]) & np.isfinite(obs_h[(s, "hap2")][ix])
        j = ix[k]
        if len(j) < 30:
            return np.nan, int(k.sum())
        a, b = obs_h[(s, "hap1")][j], obs_h[(s, "hap2")][j]
        if np.ptp(a) == 0 or np.ptp(b) == 0:
            return np.nan, int(k.sum())
        return float(spearmanr(a, b).statistic), int(k.sum())

    # ---- 点估计 ----
    rows = []
    allix = np.arange(len(chrom))
    for s in samples:
        r, n = rs_cov(s, allix)
        rf = sb(r)
        # 有效性规则：两条单倍型都没有覆盖的拷贝过半（不是「双零过半」——双零还包含
        # 有覆盖但确实关闭的拷贝，那是可测的）
        a1, a2, m = hap[f"{s}.hap1"].values, hap[f"{s}.hap2"].values, pool[s].values
        frac0 = float(np.mean((a1 == 0) & (a2 == 0)))
        c1, c2 = hcov[f"{s}.hap1"].values, hcov[f"{s}.hap2"].values
        frac_nocov = float(np.mean((c1 == 0) & (c2 == 0)))
        rho1 = spearmanr(m, a1).statistic
        rho2 = spearmanr(m, a2).statistic
        valid = not (frac_nocov > 0.5 or min(rho1, rho2) < 0.2)
        note = ""
        if frac_nocov > 0.5:
            note = "both haplotypes uncovered in >50% of copies"
        elif min(rho1, rho2) < 0.2:
            note = "merged-hap rho<0.2"
        grp = ("pluripotent" if s in PLURI else
               ("other (LCL)" if str(man.set_index("sample").loc[s, "tissue"]).startswith("Lymphoblastoid")
                else "other"))
        rows.append(dict(sample=s, group=grp, n_cov=n, frac_both0=round(frac0, 4),
                         frac_nocov=round(frac_nocov, 4),
                         rho_m_h1=round(float(rho1), 4), rho_m_h2=round(float(rho2), 4),
                         med_pa5_cov=round(float(np.nanmedian(obs_pool[s])), 3),
                         r_s_cov=None if not np.isfinite(r) else round(r, 4),
                         r_full_cov=None if not np.isfinite(rf) else round(rf, 4),
                         valid=bool(valid), note=note))
    df = pd.DataFrame(rows)

    # ---- 块 bootstrap（共同随机数） ----
    cats = pd.Categorical(chrom)
    gix = [np.where(cats.codes == j)[0] for j in range(len(cats.categories))]
    rng = np.random.default_rng(SEED)
    boots = [np.concatenate([gix[j] for j in rng.integers(0, len(gix), len(gix))])
             for _ in range(B)]
    boot = {s: {"rf": [], "rho": [], "pair": []} for s in samples}
    for b, ix in enumerate(boots):
        rfh = sb(rs_cov("H1", ix)[0])
        boot["H1"]["rf"].append(rfh)
        mh1 = obs_pool["H1"][ix]
        for s in samples:
            rf = sb(rs_cov(s, ix)[0])
            boot[s]["rf"].append(rf)
            if s == "H1":
                continue
            ms = obs_pool[s][ix]
            okp = np.isfinite(ms) & np.isfinite(mh1)
            rho = float(spearmanr(ms[okp], mh1[okp]).statistic) if okp.sum() > 30 and \
                np.ptp(ms[okp]) > 0 and np.ptp(mh1[okp]) > 0 else np.nan
            boot[s]["rho"].append(rho)
            boot[s]["pair"].append((rho, rf, rfh))
        if (b + 1) % 200 == 0:
            print(f"  bootstrap {b+1}/{B}", flush=True)

    def ci(v):
        v = np.array([x for x in v if np.isfinite(x)])
        return (np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5)) if len(v) > 50 else (np.nan, np.nan)

    out = []
    for _, r in df.iterrows():
        s = r["sample"]
        lo, hi = ci(boot[s]["rf"])
        rf = r["r_full_cov"]
        if not bool(r["valid"]) or not np.isfinite(rf):
            read = "n/a"
        elif rf >= 0.667:
            read = "borderline" if (np.isfinite(lo) and lo < 0.667) else "yes"
        else:
            read = "no"
        out.append(dict(r, r_full_cov_lo=None if not np.isfinite(lo) else round(lo, 4),
                        r_full_cov_hi=None if not np.isfinite(hi) else round(hi, 4), readable=read))
    rel = pd.DataFrame(out)
    rel.to_csv(os.path.join(OUT, "t15_reliability_cov_41.tsv"), sep="\t", index=False)

    # ---- ρ vs H1 ----
    rhos = []
    for _, r in rel.iterrows():
        s = r["sample"]
        ms, mh1 = obs_pool[s], obs_pool["H1"]
        okp = np.isfinite(ms) & np.isfinite(mh1)
        rho = float(spearmanr(ms[okp], mh1[okp]).statistic) if okp.sum() > 30 else np.nan
        rfo = r["r_full_cov"]
        rfh = float(rel.loc[rel["sample"] == "H1", "r_full_cov"].iloc[0])
        # 只在可靠度不至于把校正放大的样本上报 ρ_corr（rf < 0.2 时 1/√rf 放大 2 倍以上）
        rc = (min(rho / np.sqrt(rfo * rfh), 1.0)
              if np.isfinite(rho) and np.isfinite(rfo) and rfo >= 0.2 and bool(r["valid"])
              else np.nan)
        pr = [x for x in boot[s]["pair"] if np.all(np.isfinite(x)) and x[1] > 0 and x[2] > 0]
        if len(pr) > 50:
            arr = np.array(pr)
            rc_b = np.clip(arr[:, 0] / np.sqrt(arr[:, 1] * arr[:, 2]), -1.0, 1.0)
            rlo, rhi = np.nanpercentile(rc_b, [2.5, 97.5])
            olo, ohi = np.nanpercentile(arr[:, 0], [2.5, 97.5])
        else:
            rlo = rhi = olo = ohi = np.nan
        rhos.append(dict(sample=s, group=r["group"], readable=r["readable"],
                         rho_obs=None if not np.isfinite(rho) else round(rho, 4),
                         rho_obs_lo=None if not np.isfinite(olo) else round(olo, 4),
                         rho_obs_hi=None if not np.isfinite(ohi) else round(ohi, 4),
                         r_sub=None if not np.isfinite(rfo) else round(rfo, 4),
                         r_ref=round(rfh, 4),
                         rho_corr=None if not np.isfinite(rc) else round(rc, 4),
                         lo=None if not np.isfinite(rlo) else round(rlo, 4),
                         hi=None if not np.isfinite(rhi) else round(rhi, 4)))
    rt = pd.DataFrame(rhos).sort_values("rho_corr", ascending=False, na_position="last")
    rt.to_csv(os.path.join(OUT, "t15_rho_vs_h1_cov_41.tsv"), sep="\t", index=False)

    print(rel[["sample", "group", "n_cov", "med_pa5_cov", "r_s_cov", "r_full_cov",
               "r_full_cov_lo", "r_full_cov_hi", "readable"]].to_string(index=False))
    print("\n可读:", int((rel.readable == "yes").sum()), " 临界:", int((rel.readable == "borderline").sum()),
          " 不可读:", int((rel.readable == "no").sum()), " n/a:", int((rel.readable == "n/a").sum()))


if __name__ == "__main__":
    main()
