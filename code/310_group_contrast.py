#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
310 — 统一组间差值口径（ROUND13 P0-1）+ 逐样本效应与留一（P0-2）。

问题
----
292 的表把三个量放在一起：
  med_pluri / med_nonpluri  = 每个样本先取拷贝中位数，再对样本求均值（sample-first）
  delta                     = 先对每个拷贝取样本均值，再取拷贝中位数（copy-first）
  did                       = sample-first，但参照集只取与焦点家族**共有**的染色体
因此 raw_diff(fam) − raw_diff(ref) = DiD 不成立（L1HS 5p1k：13.22−2.75=10.47，
减 0.39 应得 10.08，表里是 10.18）。

本轮
----
1. 全部改成 sample-first：每样本取拷贝中位数 -> 组内均值 -> 相减。
2. 参照集用**全部** L1PA4–7 拷贝；bootstrap 在家族集与参照集内**各自**按染色体块重采样。
3. 强制恒等式 raw_diff(fam) − raw_diff(ref) = DiD（断言 1e-9）。
4. copy-first 的 Δ 另表保留，列名写清。
5. P0-2：11 个样本的 L1PA2 与 L1PA4–7 5p1k 中位数、差值、留一 DiD、精确置换。

输出：t13/t13_group_contrast_sample_first.tsv
      t13/t13_group_contrast_copy_first.tsv
      t13/t13_per_sample_did.tsv
      t13/t13_group_contrast.md
"""
import os
import numpy as np
import pandas as pd
import pybigtools

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
OUT = os.path.join(RES, "t13")
RMSK = "../../01_data/ref/hg38.rmsk.txt.gz"
PA_DIR = "../../01_data/fiberseq/pa"
B = 1000
SEED = 20261124
CAP = 5000
PLURI = ["H1", "H9", "WTC-11", "GM27730"]
NONPL = ["K562", "HepG2", "Hap1", "Hek293T", "MCF-7", "HG002", "GM12878"]
SAMPLES = PLURI + NONPL

PA_LOCAL = {
    "H1": "../../Age2.3.2/01_data/fiberseq/all.percent.accessible.bw",
    "H9": "../../Age2.3.2/01_data/fiberseq/h9/H9_all.percent.accessible.bw",
    **{s: os.path.join(PA_DIR, f"{s}_all.percent.accessible.bw") for s in NONPL + ["WTC-11", "GM27730"]},
}

FAMS = {
    "LTR7":    dict(pat=r"^LTR7$", minlen=350, l1=False),
    "LTR7B":   dict(pat=r"^LTR7B$", minlen=350, l1=False),
    "LTR7Y":   dict(pat=r"^LTR7Y$", minlen=350, l1=False),
    "LTR5_Hs": dict(pat=r"^LTR5_Hs$", minlen=700, l1=False),
    "L1HS":    dict(pat=r"^L1HS$", minlen=5000, l1=True),
    "SVA_D":   dict(pat=r"^SVA_D$", minlen=1000, l1=False),
    "SVA_E":   dict(pat=r"^SVA_E$", minlen=800, l1=False),
    "SVA_F":   dict(pat=r"^SVA_F$", minlen=800, l1=False),
    "AluY":    dict(pat=r"^AluY$", minlen=250, l1=False, cap=5000),
    "AluS":    dict(pat=r"^(AluSx|AluSq|AluSg|AluSp|AluSc|AluSx1)$", minlen=250, l1=False, cap=5000),
}


def load_rmsk():
    cols = ["bin", "swScore", "milliDiv", "milliDel", "milliIns", "chrom", "start", "end",
            "genoLeft", "strand", "repName", "repClass", "repFamily", "repStart", "repEnd",
            "repLeft", "id"]
    return pd.read_csv(RMSK, sep="\t", header=None, names=cols, usecols=range(17),
                       dtype={"chrom": str})


def l1_panel():
    u = pd.read_csv(os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz"), sep="\t")
    u = u[u.part == "L1PA_5p"].drop_duplicates("unit_id").reset_index(drop=True)
    sp = u.unit_id.str.rsplit("_", n=3, expand=True)
    u["el_start"] = sp[1].astype(int); u["el_end"] = sp[2].astype(int); u["grp"] = sp[3]
    u["strand"] = np.where(np.abs(u.start - u.el_start) <= 2, "+", "-")
    u.loc[np.abs(u.end - u.el_end) <= 2, "strand"] = "-"
    out = {}
    for fam, groups in (("L1PA2", ["L1PA2"]), ("L1PA3", ["L1PA3"]),
                        ("L1PA4-7", ["L1PA4", "L1PA5", "L1PA6", "L1PA7"])):
        d = u[u.grp.isin(groups)]
        out[fam] = d[["chrom", "el_start", "el_end", "strand"]].rename(
            columns={"el_start": "start", "el_end": "end"}).reset_index(drop=True)
    return out


def build():
    r = load_rmsk()
    r = r[r.chrom.str.match(r"^chr([0-9]{1,2}|X|Y)$")]
    fams = {}
    for fam, cfg in FAMS.items():
        d = r[r.repName.str.match(cfg["pat"])].copy()
        d["elen"] = d.end - d.start
        d = d[d.elen >= cfg["minlen"]]
        if "cap" in cfg and len(d) > cfg["cap"]:
            d = d.sample(cfg["cap"], random_state=SEED)
        fams[fam] = d[["chrom", "start", "end", "strand"]].reset_index(drop=True)
    fams.update(l1_panel())
    return fams


def scan_family(d, samples, cap=None):
    """与 292 完全相同的读取方式：从生物 5′ 端读 max(w)，再切 400 bp / 1 kb / 全长。"""
    res = {}
    n = len(d)
    for s in samples:
        bw = pybigtools.open(PA_LOCAL[s])
        m400 = np.full(n, np.nan); mfull = np.full(n, np.nan); m1k = np.full(n, np.nan)
        for i, (c, a, b, st) in enumerate(zip(d.chrom, d.start, d.end, d.strand)):
            a, b = int(a), int(b)
            L = b - a
            if L <= 0:
                continue
            w = L if cap is None else min(L, cap)
            if st == "+":
                rs, re_ = a, a + w
            else:
                rs, re_ = b - w, b
            v = np.asarray(bw.values(c, rs, re_, fillna=np.nan), dtype=float)
            if len(v) != w:
                continue
            f = np.isfinite(v)
            v0 = np.where(f, v, 0.0)
            k4, k1 = min(400, w), min(1000, w)
            if st == "+":
                s400, s1k = v0[:k4], v0[:k1]
            else:
                s400, s1k = v0[-k4:], v0[-k1:]
            m400[i] = s400.mean(); mfull[i] = v0.mean(); m1k[i] = s1k.mean()
        bw.close()
        res[s] = {"start400": m400, "full": mfull, "5p1k": m1k}
    return res


def blocks(chrom):
    c = pd.Categorical(np.asarray(chrom, str))
    return [np.where(c.codes == j)[0] for j in range(len(c.categories))]


def sample_medians(M, ix):
    return np.array([np.nanmedian(M[ix, j]) for j in range(M.shape[1])])


def contrasts(fam_mat, ref_mat, fam_chrom, ref_chrom, B=B, seed=SEED):
    """sample-first：每样本拷贝中位数 -> 组均值差。bootstrap 两套集合各自按块重采样。"""
    fb, rb = blocks(fam_chrom), blocks(ref_chrom)

    def one(fi, ri):
        fm = sample_medians(fam_mat, fi)
        rm = sample_medians(ref_mat, ri)
        raw = float(np.nanmean(fm[:4]) - np.nanmean(fm[4:]))
        ref = float(np.nanmean(rm[:4]) - np.nanmean(rm[4:]))
        return raw, ref, raw - ref

    raw0, ref0, did0 = one(np.arange(len(fam_chrom)), np.arange(len(ref_chrom)))
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(B):
        fi = np.concatenate([fb[j] for j in rng.integers(0, len(fb), len(fb))])
        ri = np.concatenate([rb[j] for j in rng.integers(0, len(rb), len(rb))])
        vals.append(one(fi, ri))
    v = np.array(vals)
    q = lambda x: (float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5)))
    return dict(raw_diff=raw0, ref_raw_diff=ref0, did=did0,
                raw_lo=q(v[:, 0])[0], raw_hi=q(v[:, 0])[1],
                ref_lo=q(v[:, 1])[0], ref_hi=q(v[:, 1])[1],
                did_lo=q(v[:, 2])[0], did_hi=q(v[:, 2])[1],
                med_pluri=float(np.nanmean(sample_medians(fam_mat, np.arange(len(fam_chrom)))[:4])),
                med_nonpluri=float(np.nanmean(sample_medians(fam_mat, np.arange(len(fam_chrom)))[4:])))


def copyfirst_delta(M, chrom, B=B, seed=SEED):
    """copy-first：先对每个拷贝取组均值差，再取拷贝中位数（292 的 delta 列）。"""
    bb = blocks(chrom)
    f = lambda ix: float(np.nanmedian(np.nanmean(M[ix][:, :4], axis=1) -
                                      np.nanmean(M[ix][:, 4:], axis=1)))
    d0 = f(np.arange(len(chrom)))
    rng = np.random.default_rng(seed)
    vals = [f(np.concatenate([bb[j] for j in rng.integers(0, len(bb), len(bb))]))
            for _ in range(B)]
    return d0, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    os.makedirs(OUT, exist_ok=True)
    print("建家族 ...", flush=True)
    fams = build()
    cache = {}
    for f, d in fams.items():
        cap = CAP if f in ("AluY", "AluS") else None
        cache[f] = (d, scan_family(d, SAMPLES, cap=cap))
        print(f"  {f:9s} n={len(d):5d}", flush=True)

    refd, R_ref = cache["L1PA4-7"]
    rows, cf_rows = [], []
    for f, (d, R) in cache.items():
        chrom = d.chrom.to_numpy()
        for w in ("start400", "full", "5p1k"):
            if w == "5p1k" and not (FAMS.get(f, {}).get("l1") or f.startswith("L1PA")):
                continue
            M = np.column_stack([R[s][w] for s in SAMPLES])
            RM = np.column_stack([R_ref[s][w] for s in SAMPLES])
            c = contrasts(M, RM, chrom, refd.chrom.to_numpy())
            cf = copyfirst_delta(M, chrom)
            cf_rows.append(dict(family=f, window=w, n_copy=len(d),
                                copy_first_delta=cf[0], lo=cf[1], hi=cf[2],
                                note="先对每个拷贝取组均值差，再取拷贝中位数；与 DiD 不可混用"))
            if f == "L1PA4-7":
                rows.append(dict(family=f, window=w, n_copy=len(d),
                                 group_pluri=c["med_pluri"], group_other=c["med_nonpluri"],
                                 raw_diff=c["raw_diff"], raw_lo=c["raw_lo"], raw_hi=c["raw_hi"],
                                 ref_raw_diff=np.nan, ref_lo=np.nan, ref_hi=np.nan,
                                 did=np.nan, did_lo=np.nan, did_hi=np.nan,
                                 reference=True))
            else:
                rows.append(dict(family=f, window=w, n_copy=len(d),
                                 group_pluri=c["med_pluri"], group_other=c["med_nonpluri"],
                                 raw_diff=c["raw_diff"], raw_lo=c["raw_lo"], raw_hi=c["raw_hi"],
                                 ref_raw_diff=c["ref_raw_diff"], ref_lo=c["ref_lo"], ref_hi=c["ref_hi"],
                                 did=c["did"], did_lo=c["did_lo"], did_hi=c["did_hi"],
                                 reference=False))
                # 恒等式
                assert abs(c["raw_diff"] - c["ref_raw_diff"] - c["did"]) < 1e-9, (f, w)

    G = pd.DataFrame(rows)
    C = pd.DataFrame(cf_rows)
    G.to_csv(os.path.join(OUT, "t13_group_contrast_sample_first.tsv"), sep="\t", index=False)
    C.to_csv(os.path.join(OUT, "t13_group_contrast_copy_first.tsv"), sep="\t", index=False)

    # ---------- P0-2：逐样本与留一（L1PA2 vs L1PA4–7，5p1k） ----------
    M2 = np.column_stack([cache["L1PA2"][1][s]["5p1k"] for s in SAMPLES])
    M47 = np.column_stack([R_ref[s]["5p1k"] for s in SAMPLES])
    m2 = sample_medians(M2, np.arange(M2.shape[0]))
    m47 = sample_medians(M47, np.arange(M47.shape[0]))
    diff = m2 - m47
    base = diff[:4].mean() - diff[4:].mean()
    ps = []
    for j, s in enumerate(SAMPLES):
        loo = np.delete(diff, j)
        ps.append(dict(sample=s, group="pluripotent" if s in PLURI else "other",
                       l1pa2_median=float(m2[j]), l1pa4_7_median=float(m47[j]),
                       diff=float(diff[j]),
                       did_after_loo=float(loo[:3].mean() - loo[3:].mean()) if s in PLURI
                       else np.nan,
                       role="reference (H1)" if s == "H1" else ""))
    P = pd.DataFrame(ps)
    P.to_csv(os.path.join(OUT, "t13_per_sample_did.tsv"), sep="\t", index=False)

    # 精确置换：11 取 4 / 留一后 11 取 3
    from itertools import combinations
    def exact_p(vals, k, obs):
        cnt = n = 0
        for c in combinations(range(len(vals)), k):
            n += 1
            if vals[list(c)].mean() - np.delete(vals, list(c)).mean() >= obs - 1e-12:
                cnt += 1
        return cnt, n, cnt / n
    c4, n4, p4 = exact_p(diff, 4, base)
    loo_rows = []
    for j in [i for i, s in enumerate(SAMPLES) if s in PLURI]:
        v = np.delete(diff, j)
        obs = v[:3].mean() - v[3:].mean()
        c3, n3, p3 = exact_p(v, 3, obs)
        loo_rows.append(dict(dropped=SAMPLES[j], did=float(obs), hits=c3, n_perm=n3, p_exact=p3))
    L = pd.DataFrame(loo_rows)
    L.to_csv(os.path.join(OUT, "t13_per_sample_did_loo_perm.tsv"), sep="\t", index=False)

    # ---------- 报告 ----------
    def fmt(x, d=2):
        return "NA" if not np.isfinite(x) else f"{x:+.{d}f}"
    lines = ["# 统一组间差值口径（310 / P0-1、P0-2）", "",
             f"样本 {len(SAMPLES)}（多能 4、非多能 7）。口径：每样本先取拷贝中位数，组内取均值，"
             "再相减；DiD = raw_diff(家族) − raw_diff(L1PA4–7)。参照集用全部 L1PA4–7 拷贝，"
             f"bootstrap 在家族集与参照集内各自按 23 条染色体块重采样，B={B}。", "",
             "| family | window | n | pluri | other | raw_diff [95%] | ref_raw_diff [95%] | DiD [95%] |",
             "|---|---|---|---|---|---|---|---|"]
    for _, r in G.iterrows():
        lines.append("| {family} | {window} | {n} | {p} | {o} | {rd} [{rl}, {rh}] | {rr} | {dd} [{dl}, {dh}] |".format(
            family=r.family, window=r.window, n=r.n_copy, p=fmt(r.group_pluri), o=fmt(r.group_other),
            rd=fmt(r.raw_diff), rl=fmt(r.raw_lo), rh=fmt(r.raw_hi),
            rr=("—" if r.reference else f"{fmt(r.ref_raw_diff)} [{fmt(r.ref_lo)}, {fmt(r.ref_hi)}]"),
            dd=("—" if r.reference else fmt(r.did)), dl=("—" if r.reference else fmt(r.did_lo)),
            dh=("—" if r.reference else fmt(r.did_hi))))
    g = lambda f, w: G[(G.family == f) & (G.window == w)].iloc[0]
    old = {"L1HS": 10.181, "L1PA2": 12.266, "L1PA3": 1.065}
    lines += ["", "## 与 292 旧值的差（5p1k）", "",
              "| family | 旧 DiD（共有染色体） | 新 DiD（全拷贝） | 差 |", "|---|---|---|---|"]
    for f, o in old.items():
        n = g(f, "5p1k").did
        lines.append(f"| {f} | {o:+.3f} | {n:+.3f} | {n-o:+.3f} |")
    lines += ["", "## P0-2 逐样本与留一", "",
              "| sample | group | L1PA2 中位 | L1PA4–7 中位 | 差 | 留一后 DiD |", "|---|---|---|---|---|---|"]
    for _, r in P.iterrows():
        lines.append(f"| {r['sample']} | {r.group} | {r.l1pa2_median:.2f} | {r.l1pa4_7_median:.2f} | "
                     f"{r['diff']:+.2f} | {'' if not np.isfinite(r.did_after_loo) else f'{r.did_after_loo:.2f}'} |")
    lines += ["", f"基准 DiD = {base:+.3f}；精确置换（C(11,4)=330）命中 {c4}/{n4}，p = {p4:.4f}。", "",
              "| 留一 | DiD | 精确置换命中 | 分母 | p |", "|---|---|---|---|---|"]
    for _, r in L.iterrows():
        lines.append(f"| 去掉 {r.dropped} | {r.did:.2f} | {r.hits} | {r.n_perm} | {r.p_exact:.4f} |")
    open(os.path.join(OUT, "t13_group_contrast.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("\n".join(lines[-40:]))
    print("\n[done] ->", OUT, flush=True)


if __name__ == "__main__":
    main()
