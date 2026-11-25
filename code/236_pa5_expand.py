#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
236 — 把 pa5（percent.accessible）从 5 个样本扩到 9 个，并在 L1PA2 单独重算 H1 相关。

为什么做这一步
--------------
上一轮（232）报的"逐拷贝相关"是在 **7,222 个拷贝混在一起** 上算的。混合池有两个问题：
  1. 亚家族之间的基线开放度差异很大（L1PA2 与 L1PA8 差几倍），混合池的 rho 有很大一部分
     来自"亚家族排序在两株细胞里一致"，而不是"同一个拷贝在两种细胞里一致"。
  2. GM27730（LCL）在混合池的 rho=0.657 异常高，与其它 LCL（0.24–0.38）差一截。
     这个数只出现在 c5 口径里，无法判断是生物信号还是 c5 的假象。

所以本轮
--------
A. 补 4 个样本的 percent.accessible（GM27730 / Hap1 / Hek293T / MCF-7，合计 15.4 GB），
   pa5 样本数 5 -> 9。九株里 3 株多能（H1/H9/WTC-11）、6 株分化（K562/HepG2/GM27730/
   Hap1/Hek293T/MCF-7），第一次可以做"多能 vs 分化"的完整排序。
B. **全部相关都在亚家族内部算**，不再报混合的 7,222（那一行本轮撤回，见 md）。

预先写死的判据（看结果之前就定好，不得改）
------------------------------------------
GM27730 的 L1PA2 逐拷贝 rho（pa5 口径）：
  rho <= 0.50   -> C5_ARTIFACT      混合池里 0.657 是 c5 口径的假象
  0.50 < rho < 0.80 -> AMBIGUOUS    不足以支持也不足以推翻
  rho >= 0.80   -> STATE_SPECIFICITY_DOWNGRADED
                 （一株分化细胞与 H1 的相关达到多能细胞水平，"状态特异"要么降级为
                  "细胞系特异"，要么承认逐拷贝模式在 LCL 里也能复现）

状态特异性的操作性定义（也用 pa5 的 L1PA2 子集）：
  H1 对 8 个其它样本的 rho 排名，若两株多能细胞（H9、WTC-11）正好占前两名，
  则 STATE_SPECIFIC_RANK 成立；否则不成立。这条比"某个阈值"稳，避开了量纲问题。

measure 列
----------
  percent_accessible   41 样本 compendium 的 all.percent.accessible.bw（0-100），5' 1 kb 窗均值
  c5_fire_coverage     fsp_panel_copy 的 c5（FIRE fiber 区间覆盖比例，0-1），无阈值
"""
import os, itertools
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
import pybigtools

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
PA_DIR = "../01_data/fiberseq/pa"
WIN = 1000
SEED = 20261124
PLURI = ["H1", "H9", "WTC-11"]
DIFF = ["K562", "HepG2", "GM27730", "Hap1", "Hek293T", "MCF-7"]

PA_LOCAL = {
    "H1":      "../Age2.3.2/01_data/fiberseq/all.percent.accessible.bw",
    "H9":      "../Age2.3.2/01_data/fiberseq/h9/H9_all.percent.accessible.bw",
    "WTC-11":  os.path.join(PA_DIR, "WTC-11_all.percent.accessible.bw"),
    "K562":    os.path.join(PA_DIR, "K562_all.percent.accessible.bw"),
    "HepG2":   os.path.join(PA_DIR, "HepG2_all.percent.accessible.bw"),
    "GM27730": os.path.join(PA_DIR, "GM27730_all.percent.accessible.bw"),
    "Hap1":    os.path.join(PA_DIR, "Hap1_all.percent.accessible.bw"),
    "Hek293T": os.path.join(PA_DIR, "Hek293T_all.percent.accessible.bw"),
    "MCF-7":   os.path.join(PA_DIR, "MCF-7_all.percent.accessible.bw"),
}
OK_GM = (0.50, 0.80)   # 判据边界，写死


def load_units():
    t = pd.read_csv(os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz"), sep="\t",
                    usecols=["unit_id", "part", "chrom", "start", "end"])
    p5 = t[t.part == "L1PA_5p"].drop_duplicates("unit_id").set_index("unit_id")
    cp = pd.read_csv(os.path.join(RES, "fsp_panel_copy.tsv.gz"), sep="\t",
                     usecols=["unit_id", "grp"])
    units = cp.drop_duplicates("unit_id").reset_index(drop=True)
    units["chrom"] = units.unit_id.map(p5.chrom)
    units["s5"] = units.unit_id.map(p5.start).astype(np.int64)
    units["e5"] = units.unit_id.map(p5.end).astype(np.int64)
    assert units.s5.notna().all()
    return units


def pa5_for(path, units, label):
    try:
        bw = pybigtools.open(path)
    except Exception as e:
        print(f"  {label:8s} 打不开（多半是下载未完成）：{type(e).__name__} —— 跳过，不猜", flush=True)
        return None
    ch = bw.chroms()
    out = np.full(len(units), np.nan, np.float64)
    for chrom, idx in units.groupby("chrom").indices.items():
        if chrom not in ch:
            continue
        for i in idx:
            a = int(units.s5.iat[i]); b = int(units.e5.iat[i])
            if b - a != WIN or b > ch[chrom]:
                continue
            try:
                v = np.asarray(bw.values(chrom, a, b, fillna=0.0), dtype=np.float64)
            except Exception:
                continue
            if len(v) == WIN:
                out[i] = v.mean()
    bw.close()
    print(f"  {label:8s} 覆盖 {np.isfinite(out).sum():,}/{len(units):,}  "
          f"均值 {np.nanmean(out):.3f}  非零 {np.nanmean(out > 0)*100:.1f}%", flush=True)
    return out


def _rank(v):
    """平均秩（处理并列）。"""
    order = np.argsort(v, kind="mergesort")
    out = np.empty(len(v), np.float64)
    out[order] = np.arange(1, len(v) + 1, dtype=np.float64)
    # 并列取平均
    sv = v[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            out[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return out


def _rho(x, y):
    rx = _rank(x); ry = _rank(y)
    rx -= rx.mean(); ry -= ry.mean()
    d = np.sqrt((rx @ rx) * (ry @ ry))
    return float(rx @ ry / d) if d > 0 else np.nan


def spearman_block(x, y, chrom, B=1000, seed=SEED):
    """按染色体的块 bootstrap。快实现：自己算秩，避免 pandas 每次调用的开销。"""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30:
        return np.nan, np.nan, np.nan, int(ok.sum())
    x, y, c = x[ok], y[ok], np.asarray(chrom)[ok]
    rho = _rho(x, y)
    cats = pd.Categorical(c)
    G = len(cats.categories)
    gix = [np.where(cats.codes == j)[0] for j in range(G)]
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(B):
        pick = rng.integers(0, G, G)
        ix = np.concatenate([gix[j] for j in pick])
        if len(ix) < 30:
            continue
        v = _rho(x[ix], y[ix])
        if np.isfinite(v):
            bs.append(v)
    lo, hi = np.percentile(bs, [2.5, 97.5]) if len(bs) > 20 else (np.nan, np.nan)
    return rho, float(lo), float(hi), int(ok.sum())


def main():
    units = load_units()
    CH = units.chrom.to_numpy()
    print(f"units {len(units):,}  亚家族 {units.grp.value_counts().to_dict()}", flush=True)

    # ---------- A. pa5：前 5 个复用，后 4 个新算 ----------
    print("\n=== A. pa5 逐拷贝 ===", flush=True)
    prev_path = os.path.join(RES, "t1c_pa5_by_copy.tsv.gz")
    P = pd.DataFrame({"unit_id": units.unit_id, "grp": units.grp})
    if os.path.exists(prev_path):
        old = pd.read_csv(prev_path, sep="\t")
        assert (old.unit_id.to_numpy() == units.unit_id.to_numpy()).all(), "旧 pa5 行序不一致"
        for c in old.columns:
            if c in ("unit_id", "grp"):
                continue
            P[c] = old[c].to_numpy()
        print(f"  复用 {prev_path}：{list(old.columns[2:])}", flush=True)
    for s, p in PA_LOCAL.items():
        if s in P.columns:
            continue
        if not os.path.exists(p):
            print(f"  {s}: 缺文件 {p}（跳过，不猜）", flush=True)
            continue
        v = pa5_for(p, units, s)
        if v is None:
            continue
        P[s] = v
    want = ["H1", "H9", "WTC-11", "K562", "HepG2", "GM27730", "Hap1", "Hek293T", "MCF-7"]
    cols = [c for c in want if c in P.columns]
    P = P[["unit_id", "grp"] + cols]
    P.to_csv(os.path.join(RES, "t1d_pa5_by_copy_9.tsv.gz"), sep="\t",
             index=False, compression="gzip")
    print(f"  保存 t1d_pa5_by_copy_9.tsv.gz  {P.shape}", flush=True)

    # ---------- B. c5 矩阵（41 样本） ----------
    cp = pd.read_csv(os.path.join(RES, "fsp_panel_copy.tsv.gz"), sep="\t",
                     usecols=["sample", "unit_id", "grp", "c5"])
    C = cp.pivot_table(index="unit_id", columns="sample", values="c5").reindex(units.unit_id)
    print(f"  c5 矩阵 {C.shape}", flush=True)

    # ---------- C. 亚家族内部的相关 ----------
    print("\n=== C. 亚家族内部的相关（不再报混合池）===", flush=True)
    SUBS = ["L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7", "L1PA8"]
    rows = []
    for meas, M in (("c5_fire_coverage", C), ("percent_accessible", P.set_index("unit_id"))):
        if "H1" not in M.columns:
            continue
        h1 = M["H1"].to_numpy(dtype=float)
        others = [c for c in M.columns if c not in ("H1", "grp", "measure")]
        for s in others:
            y = M[s].to_numpy(dtype=float)
            for sub in SUBS + ["L1PA13", "L1PA16", "other"]:
                m = (units.grp == sub).to_numpy()
                if m.sum() < 30:
                    continue
                rho, lo, hi, n = spearman_block(h1[m], y[m], CH[m],
                                                B=(1000 if meas == "percent_accessible" else 400))
                rows.append(dict(measure=meas, pair=f"H1_vs_{s}", other=s, subfamily=sub,
                                 n_in_subfamily=int(m.sum()), n_used=n,
                                 spearman=rho, boot_lo=lo, boot_hi=hi))
    CR = pd.DataFrame(rows)
    CR["rank_within_subfamily"] = CR.groupby(["measure", "subfamily"]).spearman.rank(ascending=False)
    CR.to_csv(os.path.join(RES, "t1d2_corr_by_subfamily.tsv"), sep="\t", index=False)
    print(f"  行 {len(CR):,}  口径 {CR.measure.unique().tolist()}", flush=True)

    # ---------- D. 判据一：GM27730 的 L1PA2 相关 ----------
    print("\n=== D. GM27730（判据一）===", flush=True)
    drow = []
    for meas in ("c5_fire_coverage", "percent_accessible"):
        s = CR[(CR.measure == meas) & (CR.other == "GM27730") & (CR.subfamily == "L1PA2")]
        if len(s) == 0:
            drow.append(dict(measure=meas, rho=np.nan, verdict="NO_DATA")); continue
        s = s.iloc[0]
        r = s.spearman
        if not np.isfinite(r):
            v = "NO_DATA"
        elif r <= OK_GM[0]:
            v = "C5_ARTIFACT"
        elif r >= OK_GM[1]:
            v = "STATE_SPECIFICITY_DOWNGRADED"
        else:
            v = "AMBIGUOUS"
        drow.append(dict(measure=meas, rho=float(r), boot_lo=s.boot_lo, boot_hi=s.boot_hi,
                         n=s.n_used, verdict=v,
                         rule=f"<={OK_GM[0]} C5_ARTIFACT; >={OK_GM[1]} DOWNGRADED"))
    DG = pd.DataFrame(drow)
    DG.to_csv(os.path.join(RES, "t1d2_gm27730_verdict.tsv"), sep="\t", index=False)
    print(DG.to_string(index=False), flush=True)

    # ---------- E. 判据二：L1PA2 子集里 H1 相关的排序 ----------
    print("\n=== E. L1PA2 子集 H1 相关的排序（判据二）===", flush=True)
    erow = []
    for meas in ("c5_fire_coverage", "percent_accessible"):
        s = CR[(CR.measure == meas) & (CR.subfamily == "L1PA2")].dropna(subset=["spearman"])
        if len(s) < 5:
            continue
        s = s.sort_values("spearman", ascending=False)
        top2 = list(s.other.head(2))
        rank_pluri = {c: int(s.index.get_loc(s.index[s.other == c][0]) + 1)
                      for c in PLURI if c in set(s.other)}
        rank_diff = {c: int(s.index.get_loc(s.index[s.other == c][0]) + 1)
                     for c in DIFF if c in set(s.other)}
        best_diff = max(rank_pluri.values()) < min(rank_diff.values()) if rank_pluri and rank_diff else False
        erow.append(dict(measure=meas, n_other=len(s), top2=",".join(top2),
                         rank_H9=rank_pluri.get("H9"), rank_WTC11=rank_pluri.get("WTC-11"),
                         rank_K562=rank_diff.get("K562"), rank_HepG2=rank_diff.get("HepG2"),
                         rank_GM27730=rank_diff.get("GM27730"),
                         rank_Hap1=rank_diff.get("Hap1"),
                         rank_Hek293T=rank_diff.get("Hek293T"),
                         rank_MCF7=rank_diff.get("MCF-7"),
                         verdict="STATE_SPECIFIC_RANK" if best_diff else "NOT_STATE_SPECIFIC_RANK",
                         top1_rho=float(s.spearman.iloc[0]), top1=top2[0] if top2 else ""))
    E = pd.DataFrame(erow)
    E.to_csv(os.path.join(RES, "t1d2_state_rank.tsv"), sep="\t", index=False)
    print(E.to_string(index=False), flush=True)

    # ---------- F. c5 与 pa5 在 L1PA2 上的对照 ----------
    print("\n=== F. c5 vs pa5（L1PA2，同一样本）===", flush=True)
    frow = []
    for s in [c for c in P.columns if c not in ("unit_id", "grp", "H1")]:
        a = CR[(CR.measure == "c5_fire_coverage") & (CR.other == s) & (CR.subfamily == "L1PA2")]
        b = CR[(CR.measure == "percent_accessible") & (CR.other == s) & (CR.subfamily == "L1PA2")]
        if len(a) == 0 or len(b) == 0:
            continue
        frow.append(dict(sample=s, rho_c5=float(a.spearman.iloc[0]), rho_pa5=float(b.spearman.iloc[0]),
                         diff_pa5_minus_c5=float(b.spearman.iloc[0] - a.spearman.iloc[0]),
                         n=a.n_used.iloc[0]))
    F = pd.DataFrame(frow).sort_values("rho_pa5", ascending=False)
    F.to_csv(os.path.join(RES, "t1d2_c5_vs_pa5_l1pa2.tsv"), sep="\t", index=False)
    print(F.to_string(index=False), flush=True)

    print("\ndone", flush=True)


if __name__ == "__main__":
    main()
