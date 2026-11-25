#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
303 — 单倍型轨道的「双零」问题（任务书 ROUND12 §0.1）。

问题
----
r_s 是 hap1 与 hap2 的 Spearman。两条轨道在**同一批拷贝**上同时为 0（该拷贝没有
定相读长、轨道被 0 填充）时，Spearman 会把它们记成「位次一致」，从而**抬高** r_s。
这些拷贝不是「真的关闭」：它们的合并轨道 pa5 与其余拷贝几乎相同
（HepG2 4.06 对 4.07，Panc1 1.77 对 1.76，Caco-2 3.83 对 4.10）。

Hap1 最极端：90.3% 的拷贝双零，ρ(合并, hap1) = −0.005，r_s = 0.79 完全没有信息。

本轮做法
--------
1. 对每个样本（L1PA2）与每个家族 × 样本，报：
   frac_both0、r_s_all、r_s_nz（去掉双零拷贝）、ρ(合并, hap1/hap2)、
   双零拷贝与其余拷贝的合并 pa5 中位。
2. 主口径改为 **r_s_nz**，r_s_all 作敏感性。
   有效性规则（先写死）：frac_both0 > 0.5 或 min(ρ(合并,hap1), ρ(合并,hap2)) < 0.2
   -> 该格 NA，note 写 `phase_missing`。
3. 用 r_s_nz 重算 r_full、校正 ρ（块 bootstrap：**每个重复内**重算 r_s_nz、
   r_H1_nz、ρ_obs 再校正，B=1000，23 条染色体块）、噪声带（σ 用新 r_full）、
   readable_0p667 的计数与分母。
4. 家族的合并轨道：L1PA2 用 t9b 的 pa 列；L1HS/LTR7/AluS 用**本地**合并 bigWig
   （H1/H9 在 Age2.3.2，其余在 01_data/fiberseq/pa），与 297 的 400 bp 窗口对齐。

判据（先写死）
--------------
R1  Spearman(pa5 中位, r_s_nz) > Spearman(pa5 中位, r_s_all)（去双零后趋势更干净）。
R2  解 r_full_nz = 0.667 的水平阈值；低于该水平的样本逐拷贝排序记为不可读。
R3  观测 ρ 落在模拟带内的样本比例；带内则「状态特异排序」不成立。
G1  多能组校正 ρ 的**区间下界最小值** > 非多能组**区间上界最大值**（不含 Hap1）。
"""
import os
import sys
import numpy as np
import pandas as pd
import pybigtools
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T9B = os.path.join(RES, "t9b")
T9H = os.path.join(RES, "t9h")
T10 = os.path.join(RES, "t10")
RMSK = "../../01_data/ref/hg38.rmsk.txt.gz"
PA_DIR = "../../01_data/fiberseq/pa"
B = 1000
B_SIM = 200
SEED = 20261124
PLURI = ["H1", "H9", "WTC-11", "GM27730"]
PA_FIXED = os.path.join(T9B, "t9b_pa5_L1PA2_41_fixed.tsv.gz")
HAP_FIXED = os.path.join(T9B, "t9b_pa5_hap_L1PA2_fixed.tsv.gz")
T1D = os.path.join(RES, "t1d_pa5_by_copy_9.tsv.gz")
FAM_PC = os.path.join(T10, "t10_family_hap_by_copy.tsv.gz")
FAM_TAB = os.path.join(T10, "t10_readability_extend.tsv")

# 家族窗口：与 297 一致（生物 5′ 端起 400 bp）
FAMS = {
    "L1HS": dict(pat=r"^L1HS$", minlen=5000),
    "LTR7": dict(pat=r"^LTR7$", minlen=350),
    "AluS": dict(pat=r"^(AluSx|AluSq|AluSg|AluSp|AluSc|AluSx1)$", minlen=250, cap=5000),
}
NMAX = 700
SUB = 400
WIN = 1000

PA_LOCAL = {
    "H1": "../../Age2.3.2/01_data/fiberseq/all.percent.accessible.bw",
    "H9": "../../Age2.3.2/01_data/fiberseq/h9/H9_all.percent.accessible.bw",
    **{s: os.path.join(PA_DIR, f"{s}_all.percent.accessible.bw")
       for s in ["WTC-11", "GM27730", "K562", "HepG2", "Hap1", "Hek293T", "MCF-7"]},
}


def sb(r):
    return 2.0 * r / (1.0 + r) if np.isfinite(r) and r > -1 else np.nan


def sp(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def block_index(chrom, rng):
    cats = pd.Categorical(chrom)
    gix = [np.where(cats.codes == j)[0] for j in range(len(cats.categories))]
    pick = rng.integers(0, len(gix), len(gix))
    return np.concatenate([gix[j] for j in pick])


def hap_stats(a1, a2, merged):
    """一个格子的双零诊断。"""
    a1 = np.asarray(a1, float); a2 = np.asarray(a2, float)
    m = np.asarray(merged, float)
    ok = np.isfinite(a1) & np.isfinite(a2)
    b0 = ok & (a1 == 0) & (a2 == 0)
    nz = ok & ~b0
    r_all = sp(a1[ok], a2[ok])
    r_nz = sp(a1[nz], a2[nz])
    return dict(
        n=int(ok.sum()), n_both0=int(b0.sum()), frac_both0=float(b0.sum() / max(ok.sum(), 1)),
        r_s_all=r_all, r_s_nz=r_nz, n_nz=int(nz.sum()),
        rho_merged_hap1=sp(m[ok], a1[ok]), rho_merged_hap2=sp(m[ok], a2[ok]),
        med_pa5_both0=float(np.nanmedian(m[b0])) if b0.sum() >= 5 else np.nan,
        med_pa5_rest=float(np.nanmedian(m[nz])) if nz.sum() >= 5 else np.nan,
        med_pa5_all=float(np.nanmedian(m[ok])) if ok.sum() else np.nan,
        cov_both0=float(np.mean(b0)) if ok.sum() else np.nan)


def phase_valid(s):
    """有效性规则（先写死）：frac_both0 > 0.5 或任一 hap 与合并轨道 ρ < 0.2 -> 无效。"""
    # 先判双零与 ρ：它们比 r_s 是否可算更能说明「这个格子没有信息」
    if s["frac_both0"] > 0.5:
        return False, "phase_missing: frac_both0>0.5"
    lo = min(np.nan_to_num(s["rho_merged_hap1"], nan=-1),
             np.nan_to_num(s["rho_merged_hap2"], nan=-1))
    if lo < 0.2:
        return False, "phase_missing: rho(merged,hap)<0.2"
    if not np.isfinite(s["r_s_nz"]):
        return False, "r_s_nan: 去双零后值退化"
    # 信度不可能是负数（它是方差比）；负估计 = 两条轨道排序相反 = 没有信息
    if s["r_s_nz"] <= 0:
        return False, "phase_conflict: r_s_nz<=0"
    return True, ""


def load_fams():
    cols = ["bin", "swScore", "milliDiv", "milliDel", "milliIns", "chrom", "start", "end",
            "genoLeft", "strand", "repName", "repClass", "repFamily", "repStart", "repEnd",
            "repLeft", "id"]
    r = pd.read_csv(RMSK, sep="\t", header=None, names=cols, usecols=range(17),
                    dtype={"chrom": str})
    r = r[r.chrom.str.match(r"^chr([0-9]{1,2}|X|Y)$")]
    out = {}
    for fam, cfg in FAMS.items():
        d = r[r.repName.str.match(cfg["pat"])].copy()
        d = d[(d.end - d.start) >= cfg["minlen"]]
        if len(d) > NMAX:
            d = d.sample(NMAX, random_state=SEED)
        out[fam] = d[["chrom", "start", "end", "strand"]].reset_index(drop=True)
    return out


def fam_windows(d):
    s = d.start.to_numpy(np.int64); e = d.end.to_numpy(np.int64)
    plus = (d.strand.to_numpy() == "+")
    rs = np.where(plus, s, np.maximum(e - WIN, 0))
    re_ = np.where(plus, s + WIN, e)
    return rs, re_, plus


def fam_merged(fam_df, rs, plus, bw):
    """合并轨道在 400 bp 窗口上的均值（与 297 的 pa400 同窗口）。"""
    n = len(fam_df)
    out = np.full(n, np.nan)
    for i, r in enumerate(fam_df.itertuples()):
        try:
            v = np.asarray(bw.values(r.chrom, int(rs[i]), int(rs[i]) + WIN, fillna=np.nan),
                           dtype=float)
        except Exception:
            continue
        f = np.isfinite(v)
        v0 = np.where(f, v, 0.0)
        seg = v0[:SUB] if plus[i] else v0[-SUB:]
        out[i] = float(seg.mean())
    return out


def main():
    os.makedirs(T9H, exist_ok=True)
    pa = pd.read_csv(PA_FIXED, sep="\t")
    hap = pd.read_csv(HAP_FIXED, sep="\t")
    units = pd.read_csv(os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz"), sep="\t")
    units = units[units.part == "L1PA_5p"].drop_duplicates("unit_id")

    # 断言（§0.5）
    t = pd.read_csv(T1D, sep="\t")
    t = t[t.grp == "L1PA2"].set_index("unit_id")
    jj = pa.set_index("unit_id").join(t[["H1"]].rename(columns={"H1": "H1_ref"}), how="inner")
    okj = np.isfinite(jj.H1) & np.isfinite(jj.H1_ref)
    rho_ck = float(spearmanr(jj.H1[okj], jj.H1_ref[okj]).statistic)
    md_ck = float(np.max(np.abs(jj.H1[okj] - jj.H1_ref[okj])))
    print(f"[assert 303] n={int(okj.sum())} rho={rho_ck:.9f} maxdiff={md_ck:.3e}", flush=True)
    assert abs(rho_ck - 1) < 1e-9 and md_ck < 1e-9, "unit_id 连接错位，中止"

    pa = pa.set_index("unit_id")
    units = units.set_index("unit_id").reindex(pa.index).reset_index()
    assert units.chrom.notna().all(), "units 缺 unit_id"
    hap = hap.set_index("unit_id").reindex(pa.index)
    chrom = units.chrom.to_numpy()
    cov = pd.read_csv(os.path.join(T9B, "t9b_coverage_verdict.tsv"), sep="\t").set_index("sample")

    # ---------- 1. L1PA2：每样本双零诊断 ----------
    cells = []
    for s in pa.columns:
        k1, k2 = f"{s}.hap1", f"{s}.hap2"
        if k1 not in hap.columns:
            continue
        st = hap_stats(hap[k1].to_numpy(float), hap[k2].to_numpy(float), pa[s].to_numpy(float))
        st.update(sample=s, family="L1PA2", is_pluri=bool(s in PLURI),
                  pa5_median=float(np.nanmedian(pa[s])),
                  n_peaks=float(cov.n_peaks.get(s, np.nan)),
                  det_pa5=float(cov.det_pa5.get(s, np.nan)),
                  frac_ge5=float(np.nanmean(pa[s] >= 5)),
                  frac_gt0=float(np.nanmean(pa[s] > 0)))
        v, note = phase_valid(st)
        st.update(phase_valid=v, note=note)
        cells.append(st)
    L1 = pd.DataFrame(cells).sort_values("pa5_median", ascending=False)
    print(L1[["sample", "frac_both0", "r_s_all", "r_s_nz", "rho_merged_hap1",
              "rho_merged_hap2", "phase_valid"]].round(3).to_string(index=False), flush=True)

    # ---------- 2. 家族 × 样本：合并轨道取本地 bigWig ----------
    fams = load_fams()
    pc = pd.read_csv(FAM_PC, sep="\t")
    fam_rows = []
    for fam, fam_df in fams.items():
        rs, re_, plus = fam_windows(fam_df)
        for s in sorted(pc[pc.family == fam]["sample"].unique()):
            sub = pc[(pc.family == fam) & (pc["sample"] == s)]   # pc.sample 是 DataFrame 方法，不能用属性取列
            assert len(sub) == len(fam_df), f"{fam} {s}: 行数不符 {len(sub)} vs {len(fam_df)}"
            if s not in PA_LOCAL or not os.path.exists(PA_LOCAL[s]):
                fam_rows.append(dict(family=fam, sample=s, note="no_merged_track"))
                continue
            bw = pybigtools.open(PA_LOCAL[s])
            merged = fam_merged(fam_df, rs, plus, bw)
            bw.close()
            st = hap_stats(sub.pa400_hap1.to_numpy(float), sub.pa400_hap2.to_numpy(float), merged)
            st.update(family=fam, sample=s, is_pluri=bool(s in PLURI),
                      pa5_median=float(np.nanmedian(merged)))
            v, note = phase_valid(st)
            st.update(phase_valid=v, note=note)
            fam_rows.append(st)
            print(f"  {fam:5s} {s:8s} both0={st['frac_both0']:.3f} "
                  f"r_all={st['r_s_all']:+.3f} r_nz={st['r_s_nz']:+.3f} "
                  f"ρ(m,h1)={st['rho_merged_hap1']:+.3f}", flush=True)
    F = pd.DataFrame(fam_rows)

    ALL = pd.concat([L1, F], ignore_index=True, sort=False)
    ALL.to_csv(os.path.join(T9H, "t10_hap_doublezero.tsv"), sep="\t", index=False)

    # ---------- 3. 用 r_s_nz 重算 r_full / 校正 ρ ----------
    h1 = L1[L1["sample"] == "H1"].iloc[0]
    paH1 = pa["H1"].to_numpy(float)
    h1a = hap["H1.hap1"].to_numpy(float)
    h1b = hap["H1.hap2"].to_numpy(float)
    rng = np.random.default_rng(SEED)
    boots = [block_index(chrom, rng) for _ in range(B)]
    print(f"块 bootstrap：{len(set(chrom))} 条染色体，B={B}", flush=True)

    corr, sim = [], []
    for _, x in L1.iterrows():
        s = x["sample"]
        if s == "H1" or not x.phase_valid:
            continue
        a1 = hap[f"{s}.hap1"].to_numpy(float)
        a2 = hap[f"{s}.hap2"].to_numpy(float)
        tgt = pa[s].to_numpy(float)
        rho_obs = sp(paH1, tgt)
        r_full_nz = sb(x.r_s_nz)
        r_full_all = sb(x.r_s_all)
        den = np.sqrt(max(r_full_nz, 1e-6) * max(sb(h1.r_s_nz), 1e-6))
        den_all = np.sqrt(max(r_full_all, 1e-6) * max(sb(h1.r_s_all), 1e-6))

        rb, cb, cba, rfb = [], [], [], []
        for ix in boots:
            v = sp(paH1[ix], tgt[ix])
            b1ix, b2ix = a1[ix], a2[ix]
            nzm = ~((b1ix == 0) & (b2ix == 0))
            rs_b = sp(b1ix[nzm], b2ix[nzm])
            h1z = ~((h1a[ix] == 0) & (h1b[ix] == 0))
            rh_b = sp(h1a[ix][h1z], h1b[ix][h1z])
            d = np.sqrt(max(sb(rs_b), 1e-6) * max(sb(rh_b), 1e-6))
            rb.append(v)
            rfb.append(sb(rs_b))
            cb.append(v / d if np.isfinite(v) and np.isfinite(d) and d > 0 else np.nan)
            rs_a = sp(b1ix, b2ix)
            rh_a = sp(h1a[ix], h1b[ix])
            da = np.sqrt(max(sb(rs_a), 1e-6) * max(sb(rh_a), 1e-6))
            cba.append(v / da if np.isfinite(v) and np.isfinite(da) and da > 0 else np.nan)
        q = lambda z: (float(np.percentile(z, 2.5)), float(np.percentile(z, 97.5))) \
            if len(z) else (np.nan, np.nan)
        rb = np.array([v for v in rb if np.isfinite(v)], float)
        cb = np.array([v for v in cb if np.isfinite(v)], float)
        cba = np.array([v for v in cba if np.isfinite(v)], float)
        rfb = np.array([v for v in rfb if np.isfinite(v)], float)
        lo, hi = q(rb); clo, chi = q(cb); alo, ahi = q(cba)
        flo, fhi = q(rfb)

        corr.append(dict(
            sample=s, rho_obs=rho_obs, lo=lo, hi=hi,
            r_s=x.r_s_nz, r_s_all=x.r_s_all, r_H1=h1.r_s_nz,
            r_full=r_full_nz, r_full_all=r_full_all, r_full_H1=sb(h1.r_s_nz),
            r_full_lo=flo, r_full_hi=fhi,
            rho_corr_sb=rho_obs / den if den > 0 else np.nan, corr_lo=clo, corr_hi=chi,
            rho_corr_sb_all=rho_obs / den_all if den_all > 0 else np.nan,
            corr_all_lo=alo, corr_all_hi=ahi,
            readable_0p5=bool(x.r_s_nz >= 0.5 and h1.r_s_nz >= 0.5),
            readable_0p667=bool(r_full_nz >= 0.667 and sb(h1.r_s_nz) >= 0.667),
            is_pluri=bool(s in PLURI)))

        # 噪声带（σ 用新的 r_full）
        h9 = pa["H9"].to_numpy(float)
        m = np.isfinite(h9) & np.isfinite(tgt)
        qtl = (pd.Series(h9[m]).rank().to_numpy() - 0.5) / m.sum()
        mapped = np.quantile(tgt[m], qtl)
        var = float(np.var(mapped)) or 1.0
        sigma = np.sqrt(max(var * (1.0 / max(r_full_nz, 0.05) - 1.0), 1e-9))
        varH1 = float(np.var(paH1[m])) or 1.0
        sigmaH1 = np.sqrt(max(varH1 * (1.0 / max(sb(h1.r_s_nz), 0.05) - 1.0), 1e-9))
        g = np.random.default_rng(SEED + 7)
        s1, s2 = [], []
        for _ in range(B_SIM):
            noisy = mapped + g.normal(0, sigma, len(mapped))
            s1.append(sp(paH1[m], noisy))
            s2.append(sp(paH1[m] + g.normal(0, sigmaH1, m.sum()), noisy))
        s1 = np.array([v for v in s1 if np.isfinite(v)])
        s2 = np.array([v for v in s2 if np.isfinite(v)])
        obs_m = sp(paH1[m], tgt[m])
        b1 = (np.percentile(s1, 2.5), np.percentile(s1, 97.5))
        b2 = (np.percentile(s2, 2.5), np.percentile(s2, 97.5))
        sim.append(dict(sample=s, obs_rho=obs_m,
                        band_lo=float(b1[0]), band_hi=float(b1[1]),
                        sim_med=float(np.median(s1)),
                        in_band=bool(b1[0] <= obs_m <= b1[1]),
                        direction=("inside" if b1[0] <= obs_m <= b1[1]
                                   else ("below" if obs_m < b1[0] else "above")),
                        sigma=float(sigma), r_full=float(r_full_nz),
                        band_lo_2side=float(b2[0]), band_hi_2side=float(b2[1]),
                        in_band_2side=bool(b2[0] <= obs_m <= b2[1]),
                        sigma_H1=float(sigmaH1)))
        print(f"  {s:8s} ρ={obs_m:+.3f} 带 [{b1[0]:+.3f},{b1[1]:+.3f}] "
              f"校正 {corr[-1]['rho_corr_sb']:+.3f} "
              f"[{corr[-1]['corr_lo']:+.3f},{corr[-1]['corr_hi']:+.3f}]", flush=True)

    C = pd.DataFrame(corr).sort_values("rho_obs", ascending=False)
    S = pd.DataFrame(sim)
    C.to_csv(os.path.join(T9H, "t9h_attenuation.tsv"), sep="\t", index=False)
    S.to_csv(os.path.join(T9H, "t9h_simulation.tsv"), sep="\t", index=False)

    # ---------- 4. 图脚本要的表 ----------
    H2 = L1.drop(columns=["med_pa5_all"]).copy()   # pa5_median 已有，避免重名
    H2["n_copy"] = len(pa)
    fig = H2.merge(C[["sample", "r_full", "r_full_all", "r_full_lo", "r_full_hi",
                      "rho_obs", "lo", "hi",
                      "rho_corr_sb", "corr_lo", "corr_hi",
                      "rho_corr_sb_all", "corr_all_lo", "corr_all_hi",
                      "readable_0p5", "readable_0p667"]],
                   on="sample", how="left")
    # H1 是参照（不在 C 里），信度本身仍要报
    hv0 = fig["sample"] == "H1"
    fig.loc[hv0, "r_full"] = sb(h1.r_s_nz)
    fig.loc[hv0, "r_full_all"] = sb(h1.r_s_all)
    # 旧口径（r_s_all）下的可读性，用于报「谁由可读变不可读」
    fig["readable_0p667_all"] = (fig.r_full_all >= 0.667) & (sb(h1.r_s_all) >= 0.667)
    fig["readable_0p5_all"] = (fig.r_s_all >= 0.5) & (h1.r_s_all >= 0.5)
    fig["flip_readable"] = (fig.readable_0p667_all.fillna(False).astype(bool)
                            & ~fig.readable_0p667.fillna(False).astype(bool))
    fig = fig.rename(columns={"lo": "rho_obs_lo", "hi": "rho_obs_hi"})
    fig["r_s"] = fig["r_s_nz"]        # 主口径
    fig["note"] = fig["note"].fillna("")
    fig.loc[fig["sample"] == "H1", "note"] = "reference"
    # H1 是参照，不参与校正（rho_corr 为 NaN），但它的信度本身可判读
    hv = fig["sample"] == "H1"
    fig.loc[hv, "readable_0p5"] = bool(h1.r_s_nz >= 0.5)
    fig.loc[hv, "readable_0p667"] = bool(sb(h1.r_s_nz) >= 0.667)
    # 无效格子（phase_missing / phase_conflict）一律不可读，且主口径 r_s 置 NA；
    # r_s_nz / r_s_all 保留原始估计，便于核对为什么被判无效
    inv = fig["phase_valid"] == False
    fig.loc[inv, ["readable_0p5", "readable_0p667"]] = False
    fig.loc[inv, "r_s"] = np.nan
    # flip 必须在 H1 与无效格子的标记**之后**算，否则 H1 会被当成「原来不可读」
    fig["flip_readable"] = (fig.readable_0p667_all.fillna(False).astype(bool)
                            & ~fig.readable_0p667.fillna(False).astype(bool))
    allpa = pa.median().rename("pmed").reset_index().rename(columns={"index": "sample"})
    extra = allpa[~allpa["sample"].isin(fig["sample"])].copy()
    extra["is_pluri"] = extra["sample"].isin(PLURI)
    extra["n_copy"] = len(pa)
    extra["n_peaks"] = [float(cov.n_peaks.get(s, np.nan)) for s in extra["sample"]]
    extra["det_pa5"] = [float(cov.det_pa5.get(s, np.nan)) for s in extra["sample"]]
    extra["note"] = "no_hap_split"
    extra = extra.rename(columns={"pmed": "pa5_median"})
    cols = ["sample", "is_pluri", "n_copy", "n_peaks", "det_pa5", "frac_ge5",
            "pa5_median", "r_s", "r_s_all", "r_s_nz", "r_full", "r_full_all",
            "rho_obs", "rho_obs_lo", "rho_obs_hi", "rho_corr_sb", "corr_lo", "corr_hi",
            "rho_corr_sb_all", "corr_all_lo", "corr_all_hi",
            "r_full_lo", "r_full_hi", "readable_0p667_all", "readable_0p5_all", "flip_readable",
            "frac_both0", "rho_merged_hap1", "rho_merged_hap2",
            "med_pa5_both0", "med_pa5_rest", "phase_valid",
            "readable_0p5", "readable_0p667", "note"]
    for c in cols:
        if c not in fig.columns:
            fig[c] = np.nan
        if c not in extra.columns:
            extra[c] = np.nan
    fig = pd.concat([fig.reindex(columns=cols), extra.reindex(columns=cols)],
                    ignore_index=True).sort_values("pa5_median", ascending=False)
    fig.to_csv(os.path.join(RES, "t10_reliability_fig.tsv"), sep="\t", index=False)

    S["pa5_median"] = [float(pa[s].median()) for s in S["sample"]]
    S["is_pluri"] = S["sample"].isin(PLURI)
    S.to_csv(os.path.join(RES, "t10_sim_band.tsv"), sep="\t", index=False)

    # 水平 vs 拷贝特异
    keep = fig[fig.r_full >= 0.667].copy()
    keep = keep[keep["sample"] != "H1"].sort_values("rho_corr_sb", ascending=False)
    keep[["sample", "is_pluri", "pa5_median", "r_full", "rho_obs",
          "rho_obs_lo", "rho_obs_hi", "rho_corr_sb", "corr_lo", "corr_hi"]].to_csv(
        os.path.join(RES, "t10_level_vs_copy.tsv"), sep="\t", index=False)

    # ---------- 5. 家族表（给 304 并进图谱） ----------
    F2 = F.copy()
    F2["r_full"] = [sb(v) for v in F2.r_s_nz]
    F2["r_full_all"] = [sb(v) for v in F2.r_s_all]
    F2["r_s"] = F2["r_s_nz"]
    F2.to_csv(os.path.join(T10, "t10_readability_nz.tsv"), sep="\t", index=False)

    # ---------- 6. 判据 ----------
    Lv = L1[L1.phase_valid]
    r1_all = float(spearmanr(Lv.pa5_median, Lv.r_s_all).statistic)
    r1_nz = float(spearmanr(Lv.pa5_median, Lv.r_s_nz).statistic)
    R1 = bool(r1_nz > r1_all)

    def interp(xs, ys, tgt):
        o = np.argsort(xs)
        for i in range(len(o) - 1):
            a, b = o[i], o[i + 1]
            if (ys[a] - tgt) * (ys[b] - tgt) <= 0 and ys[a] != ys[b]:
                return float(xs[a] + (tgt - ys[a]) * (xs[b] - xs[a]) / (ys[b] - ys[a]))
        return np.nan

    rf = np.array([sb(v) for v in Lv.r_s_nz])
    thr_full = interp(Lv.pa5_median.to_numpy(), rf, 0.667)
    thr_s = interp(Lv.pa5_median.to_numpy(), Lv.r_s_nz.to_numpy(), 0.50)

    n_out = int((~S.in_band).sum())
    pl = C[C.is_pluri]; non = C[~C.is_pluri]
    g1_lo = float(pl.corr_lo.min()); g1_hi = float(non.corr_hi.max())
    G1 = bool(g1_lo > g1_hi)

    # ---------- 7. 报告 ----------
    n_valid = int(L1.phase_valid.sum())
    n_read = int(fig.readable_0p667.fillna(False).sum())
    L = ["# 可靠度与衰减校正（303：单倍型「双零」修正）", "",
         f"L1PA2 的 {len(pa)} 个拷贝；{len(L1)} 个样本能做 hap 伪重复，"
         f"其中 **{n_valid}** 个通过有效性规则（`phase_valid`），"
         f"{len(L1)-n_valid} 个判为 `phase_missing`。", "",
         "## 0. 本轮改了什么", "",
         "r_s 是 hap1 与 hap2 的 Spearman。两条轨道在**同一批拷贝**上同时为 0 时"
         "（该拷贝没有定相读长、轨道被 0 填充），Spearman 会把它们记成「位次一致」，"
         "从而抬高 r_s。这些拷贝不是真的关闭——合并轨道的 pa5 与其余拷贝几乎相同。", "",
         "主口径改为 **r_s_nz**（去掉双零拷贝），r_s_all 作敏感性。"
         "有效性规则先写死：frac_both0 > 0.5 或 min(ρ(合并,hap1), ρ(合并,hap2)) < 0.2 "
         "-> 该格 NA 并记 `phase_missing`。", "",
         "## 1. 双零诊断（每样本）", "",
         "| 样本 | 双零比例 | 双零的 pa5 中位 | 其余的 pa5 中位 | r_s_all | r_s_nz | ρ(合并,hap1) | ρ(合并,hap2) | 有效 |",
         "|---|---|---|---|---|---|---|---|---|"]
    for _, x in L1.iterrows():
        L.append(f"| {x['sample']} | {x.frac_both0:.3f} | {x.med_pa5_both0:.2f} | "
                 f"{x.med_pa5_rest:.2f} | {x.r_s_all:.3f} | {x.r_s_nz:.3f} | "
                 f"{x.rho_merged_hap1:+.3f} | {x.rho_merged_hap2:+.3f} | "
                 f"{'是' if x.phase_valid else '否'} |")
    L += ["", "双零拷贝的合并 pa5 与其余拷贝几乎相同（差值多在 0.1 pp 以内），"
          "只有 H1/WTC-11/GM27730/H9 这类高水平样本的差值大一些（双零拷贝略低，"
          "但仍不是 0）——所以「双零 = 关闭」不成立，它是「该拷贝没有定相读长」。", "",
          "### 1.1 家族 × 样本", "",
          "| 家族 | 样本 | 双零比例 | r_s_all | r_s_nz | ρ(合并,hap1) | 有效 |",
          "|---|---|---|---|---|---|---|"]
    for _, x in F.iterrows():
        if not np.isfinite(x.get("r_s_nz", np.nan)):
            continue
        L.append(f"| {x.family} | {x['sample']} | {x.frac_both0:.3f} | {x.r_s_all:.3f} | "
                 f"{x.r_s_nz:.3f} | {x.rho_merged_hap1:+.3f} | "
                 f"{'是' if x.phase_valid else '否'} |")
    L += ["", "家族表只有 9 个样本有 hap 轨道（297 的取数范围），L1PA2 有 15 个；"
          "分母不同，不能混在一起报。", "",
          "## 2. 校正后 ρ（对 H1，r_s_nz 口径）", "",
          "| 样本 | r_full | ρ_obs | 95%（ρ_obs） | 校正 ρ | 95%（校正） | 旧口径（r_s_all） |",
          "|---|---|---|---|---|---|---|"]
    for _, x in C.iterrows():
        L.append(f"| {x['sample']} | {x.r_full:.3f} | {x.rho_obs:+.3f} | "
                 f"[{x.lo:+.3f}, {x.hi:+.3f}] | {x.rho_corr_sb:+.3f} | "
                 f"[{x.corr_lo:+.3f}, {x.corr_hi:+.3f}] | {x.rho_corr_sb_all:+.3f} |")
    L += ["",
          f"- 多能样本校正 ρ：{min(pl.rho_corr_sb):.3f}–{max(pl.rho_corr_sb):.3f}（n={len(pl)}）",
          f"- 非多能样本校正 ρ：{min(non.rho_corr_sb):.3f}–{max(non.rho_corr_sb):.3f}（n={len(non)}）",
          "",
          "### Go/No-go 0.1：区间口径（不含 Hap1）", "",
          f"- 多能组**区间下界的最小值** = **{g1_lo:.3f}**"
          f"（{pl.loc[pl.corr_lo.idxmin(),'sample']}）",
          f"- 非多能组**区间上界的最大值** = **{g1_hi:.3f}**"
          f"（{non.loc[non.corr_hi.idxmax(),'sample']}）",
          f"- 判据 G1（下界最小值 > 上界最大值）-> **{'PASS' if G1 else 'FAIL'}**，"
          f"余量 {g1_lo - g1_hi:+.3f}", "",
          "## 3. 可靠度由动态范围决定", "",
          f"- Spearman(pa5 中位, r_s_nz) = **{r1_nz:+.3f}**；"
          f"同一批样本用 r_s_all = **{r1_all:+.3f}** -> R1 "
          f"**{'PASS' if R1 else 'FAIL'}**（去双零后趋势更干净）",
          f"- r_s_nz = 0.50 对应 pa5 中位 ≈ **{thr_s:.2f}**；r_full = 0.667 对应 ≈ "
          f"**{thr_full:.2f}**" if np.isfinite(thr_s) and np.isfinite(thr_full) else
          f"- 阈值未能在样本范围内内插（r_s {thr_s:.2f}、r_full {thr_full:.2f}）",
          "", "## 4. 模拟：同一排序 + 噪声（σ 用新的 r_full）", "",
          "| 样本 | ρ 观测 | 带（单侧噪声） | 带内 | 带（双侧噪声） | 带内 |",
          "|---|---|---|---|---|---|"]
    for _, x in S.iterrows():
        L.append(f"| {x['sample']} | {x.obs_rho:+.3f} | [{x.band_lo:+.3f}, {x.band_hi:+.3f}] | "
                 f"{'是' if x.in_band else '否'} | [{x.band_lo_2side:+.3f}, {x.band_hi_2side:+.3f}] | "
                 f"{'是' if x.in_band_2side else '否'} |")
    L += ["", f"单侧噪声：**{n_out}/{len(S)}** 带外；双侧："
          f"**{int((~S.in_band_2side).sum())}/{len(S)}**。", "",
          f"## 5. 可读计数（分母要写清）", "",
          f"- 41 个样本里 15 个能做 hap 拆分；15 个里 **{n_valid}** 个通过有效性规则。",
          f"- r_s_nz ≥ 0.50 且 H1 也 ≥ 0.50 的样本：**{int(fig.readable_0p5.fillna(False).sum())}**"
          f"（分母 {len(L1)}）。",
          f"- r_full_nz ≥ 0.667 且 H1 也 ≥ 0.667 的样本：**{n_read}**（分母 {len(L1)}）。",
          f"- 旧写法「13/15 可读」**作废**：它用的是被双零抬高的 r_s_all，"
          f"且把 Hap1 算进分母。", ""]

    flips = fig[fig.flip_readable.fillna(False).astype(bool)]
    L += ["### 5.1 谁由「可读」变「不可读」", "",
          "旧口径 r_full_all ≥ 0.667、新口径 r_full_nz < 0.667 的样本：", "",
          "| 样本 | 双零比例 | r_full_all（旧） | r_full_nz（新） | 新 r_full 的 95% |",
          "|---|---|---|---|---|"]
    for _, x in flips.iterrows():
        L.append(f"| {x['sample']} | {x.frac_both0:.3f} | {x.r_full_all:.3f} | "
                 f"{x.r_full:.3f} | [{x.r_full_lo:.3f}, {x.r_full_hi:.3f}] |")
    if len(flips):
        L += ["", f"共 **{len(flips)}** 个（{', '.join(flips['sample'])}）。"
              f"双零比例 {flips.frac_both0.min():.3f}–{flips.frac_both0.max():.3f}，"
              f"r_full 被抬高 {(flips.r_full_all-flips.r_full).min():.3f}–"
              f"{(flips.r_full_all-flips.r_full).max():.3f}，去双零后掉到阈值以下。", ""]
    hg = fig[fig["sample"] == "HG002"]
    if len(hg):
        L += [f"**HG002 不在上面的名单里**，但它是同一现象的最清楚例子：r_s_all "
              f"{hg.r_s_all.iloc[0]:.3f} -> r_s_nz {hg.r_s_nz.iloc[0]:.3f}，"
              f"r_full {hg.r_full_all.iloc[0]:.3f} -> {hg.r_full.iloc[0]:.3f}。"
              "它在旧口径下也不可读（0.667 阈值），所以不产生「翻转」，"
              "但双零抬高同样存在。", ""]

    near = fig[(fig.r_s.notna()) & fig.r_full.notna()
               & (fig.r_full.between(0.60, 0.74))]
    L += ["### 5.2 阈值附近的样本（要用区间说话）", "",
          "| 样本 | r_full_nz | 95% | 是否跨 0.667 | r_s_nz |",
          "|---|---|---|---|---|"]
    for _, x in near.iterrows():
        cross = "是" if (x.r_full_lo <= 0.667 <= x.r_full_hi) else "否"
        L.append(f"| {x['sample']} | {x.r_full:.3f} | [{x.r_full_lo:.3f}, {x.r_full_hi:.3f}] | "
                 f"{cross} | {x.r_s_nz:.3f} |")
    L += ["", "跨阈值的不写「可读」或「不可读」，写「阈值附近，区间跨 0.667」。", ""]
    L += ["## 6. 限制", "",
          "- **r_full 是下界**：hap 之间的 cis 差异被算进「不可靠」。",
          "- **校正后的 ρ 是上界**：分母是下界，商即上界。只能读「至少这么接近」，"
          "不能用来排两组的细微差别。",
          "- **「双零」本身有信息损失**：删掉这些拷贝会改变样本的组成（剩下的是"
          "「有定相读长」的拷贝，可能偏向某些基因组区域）。已在第 1 节给出比例，"
          "比例高（> 20%）的样本（Caco-2 0.322、Panc1 0.285、A549 0.229、THP-1 0.196）"
          "的 r_s_nz 只能当**下界**，不能与低双零样本直接比。",
          "- **只有 15/41 个样本能做 hap 拆分**；其余 26 个的排序可读性未知。",
          "- Hap1 全部格子判 `phase_missing`（frac_both0 0.903、ρ(合并,hap1) −0.005），"
          "不得再出现在可读性趋势或图 C(d) 里。", ""]
    open(os.path.join(T9H, "t9h_reliability.md"), "w", encoding="utf-8").write("\n".join(L))
    print(f"R1 rho_nz={r1_nz:+.3f} rho_all={r1_all:+.3f} -> {R1}", flush=True)
    print(f"G1 {g1_lo:.3f} > {g1_hi:.3f} -> {G1}（余量 {g1_lo-g1_hi:+.3f}）", flush=True)
    print(f"有效 {n_valid}/{len(L1)}；readable_0p667 {n_read}", flush=True)


if __name__ == "__main__":
    main()
