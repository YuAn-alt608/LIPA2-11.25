#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
297 — 可读性图谱扩展到 L1HS / LTR7 / AluS（任务书 ROUND11 §0.4）。

为什么要做
----------
`t10_readability_map.tsv` 原来只有 L1PA2 一行有 r_s，图 C(d) 是空的。
图 C 要说明的是「可读性由水平决定，与家族无关」，那就必须有**别的家族**的点。

做法
----
对 3 个家族 × 9 个样本 × hap1/hap2，取每个拷贝生物 5′ 端起 1 kb 的
percent.accessible 均值（分母固定窗长，未覆盖按 0，与 273/292 同口径），
再从同一个数组切出前 400 bp，得到两套窗口。**每个拷贝只读一次 1 kb**，
不下载整个 bigWig（_remotebw 的 HTTP Range）。

预算（先算再跑，任务书要求）
--------------------------
单窗实测 17.6 KB、7.3 个请求（MCF-7 hap1，40 窗抽样）。
每家族取 ≤700 条（任务书上限定 1,000，这里留余量）：
    3 家族 × 700 拷贝 × 9 样本 × 2 hap = 37,800 窗
    预估 37,800 × 17.6 KB ≈ 0.65 GB，≈ 27.6 万请求
低于任务书的 1 GB 线。脚本内设 1.0 GB 硬上限：超了就停并记录哪些格子没取。

判据（先写死）
--------------
M1  L1HS / LTR7 / AluS 的 r_s 与该样本水平的 Spearman 与 L1PA2 同号且量级相近
    （同一张图上应当落在同一条趋势线上）。
M2  若某家族在某样本的 r_s 与 L1PA2 差 > 0.15，需在报告里单独说明（家族特异的
    可读性，不能当成「水平决定」）。
M3  覆盖率 < 50% 的格子标 low_coverage，不参与趋势线拟合。
"""
import os
import sys
import time
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _remotebw import RemoteBigWig

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
RMSK = "../../01_data/ref/hg38.rmsk.txt.gz"
MAN = "../../01_data/meta/fiberseq_samples.tsv"
OUT = os.path.join(RES, "t10")
SEED = 20261124
NTHREAD = 24
NMAX = 700
BYTE_CAP = 1.0e9
WIN = 1000
SUB = 400

FAMS = {
    "L1HS": dict(pat=r"^L1HS$", minlen=5000),
    "LTR7": dict(pat=r"^LTR7$", minlen=350),
    "AluS": dict(pat=r"^(AluSx|AluSq|AluSg|AluSp|AluSc|AluSx1)$", minlen=250, cap=5000),
}
SAMPLES = ["H1", "H9", "WTC-11", "GM27730", "HepG2", "K562", "Hap1", "Hek293T", "MCF-7"]
TOTAL_BYTES = [0]
TOTAL_REQ = [0]


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
        print(f"  {fam:6s} n={len(out[fam]):4d} 中位长 "
              f"{int(np.median(out[fam].end-out[fam].start))}", flush=True)
    return out


def windows(d):
    """从生物 5′ 端起 1 kb；负链取 [end-1000, end)。"""
    s = d.start.to_numpy(np.int64); e = d.end.to_numpy(np.int64)
    plus = (d.strand.to_numpy() == "+")
    rs = np.where(plus, s, np.maximum(e - WIN, 0))
    re_ = np.where(plus, s + WIN, e)
    return rs, re_, plus


def worker(url, rows, plus, retries=2):
    for attempt in range(retries + 1):
        try:
            bw = RemoteBigWig(url)
            break
        except Exception:
            if attempt == retries:
                return [(np.nan, 0.0, np.nan, 0.0)] * len(rows)
            time.sleep(1.5 * (attempt + 1))
    out = []
    try:
        for (c, a, b), pf in zip(rows, plus):
            try:
                v = np.asarray(bw.values(c, int(a), int(b), fill=np.nan), dtype=float)
                f = np.isfinite(v)
                v0 = np.where(f, v, 0.0)
                if pf:
                    s4, c4 = v0[:SUB], f[:SUB]
                else:
                    s4, c4 = v0[-SUB:], f[-SUB:]
                out.append(((float(s4.mean()), float(np.mean(c4)),
                             float(v0.mean()), float(np.mean(f)))))
            except Exception:
                out.append(((np.nan, 0.0, np.nan, 0.0)))
    finally:
        TOTAL_BYTES[0] += bw.n_bytes
        TOTAL_REQ[0] += bw.n_requests
        bw.close()
    return out


def fetch(url, rows, plus, tag=""):
    """按 stride 分块并行取数，**再按原位置放回**。

    注意：第一版把各线程的结果按块顺序拼接，等于对行做了一次步长 = 线程数的
    置换（行 0,24,48,… 排到前面）。这正是 273 产出的 t9b 文件 unit_id 与数值
    行整体错位的原因（见 298 的审计）。这里用位置数组散回去。
    """
    k = min(NTHREAD, max(4, len(rows) // 20))
    pos = np.arange(len(rows))
    chunks = [(rows[i::k], plus[i::k], pos[i::k]) for i in range(k)]
    res = [None] * k
    t0 = time.time()
    with ThreadPoolExecutor(k) as ex:
        futs = {ex.submit(worker, url, ch, pl): i for i, (ch, pl, _) in enumerate(chunks) if len(ch)}
        for f in as_completed(futs):
            res[futs[f]] = f.result()
    out = np.full((len(rows), 4), np.nan)
    for (ch, pl, ix), part in zip(chunks, res):
        if len(ch) == 0:
            continue
        assert len(part) == len(ix), f"{tag}: 块长度不符"
        out[ix] = np.array(part, float)
    arr = out
    assert np.isfinite(arr[:, 2]).sum() > 0.9 * len(rows), f"{tag}: 有效值过少"
    print(f"  {tag:26s} {len(rows):4d} 窗 {time.time()-t0:6.1f}s  中位 {np.nanmedian(arr[:,2]):6.2f}"
          f"  |  累计 {TOTAL_BYTES[0]/1e9:.3f} GB", flush=True)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def main():
    os.makedirs(OUT, exist_ok=True)
    man = pd.read_csv(MAN, sep="\t").set_index("sample")
    fams = load_fams()

    per_copy = []
    summary = []
    for fam, fam_df in fams.items():
        rs, re_, plus = windows(fam_df)
        rows = list(zip(fam_df.chrom, rs, re_))
        base = dict(family=fam)
        for s in SAMPLES:
            if s not in man.index:
                continue
            if TOTAL_BYTES[0] > BYTE_CAP:
                summary.append(dict(family=fam, sample=s, window_def="start400",
                                    r_s=np.nan, r_full=np.nan, n=0, cov_med=np.nan,
                                    note="byte_cap_reached"))
                continue
            vals = {}
            for h in ("1", "2"):
                col = f"percent_accessible_hap{h}"
                if col not in man.columns or pd.isna(man.loc[s, col]):
                    vals[h] = None
                    continue
                vals[h] = fetch(man.loc[s, col], rows, plus, tag=f"{fam} {s} hap{h}")
            if vals.get("1") is None or vals.get("2") is None:
                summary.append(dict(family=fam, sample=s, window_def="start400", r_s=np.nan,
                                    r_full=np.nan, n=0, cov_med=np.nan, note="no_hap_track"))
                continue
            a1, c1, b1, k1 = vals["1"]
            a2, c2, b2, k2 = vals["2"]
            for wd, x, y, cv in (("start400", a1, a2, (c1, c2)),
                                 ("5p1k", b1, b2, (k1, k2))):
                ok = np.isfinite(x) & np.isfinite(y)
                r = float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() > 30 else np.nan
                cov = float(np.nanmedian(np.r_[cv[0][ok], cv[1][ok]])) if ok.sum() else np.nan
                summary.append(dict(family=fam, sample=s, window_def=wd, r_s=r,
                                    r_full=2 * r / (1 + r) if np.isfinite(r) and r > -1 else np.nan,
                                    n=int(ok.sum()), cov_med=cov,
                                    note="low_coverage" if cov < 0.5 else ""))
            for i in range(len(rows)):
                per_copy.append(dict(family=fam, sample=s,
                                     pa400_hap1=a1[i], pa1k_hap1=b1[i],
                                     pa400_hap2=a2[i], pa1k_hap2=b2[i],
                                     cov400=min(c1[i], c2[i]), cov1k=min(k1[i], k2[i])))
            print(f"    {fam} {s}: r_s(400)={summary[-2]['r_s']:+.3f} "
                  f"r_s(1k)={summary[-1]['r_s']:+.3f}  bytes={TOTAL_BYTES[0]/1e9:.3f} GB",
                  flush=True)

    S = pd.DataFrame(summary)
    S.to_csv(os.path.join(OUT, "t10_readability_extend.tsv"), sep="\t", index=False)
    pd.DataFrame(per_copy).to_csv(os.path.join(OUT, "t10_family_hap_by_copy.tsv.gz"),
                                  sep="\t", index=False, compression="gzip")
    print(f"\n总字节 {TOTAL_BYTES[0]/1e9:.3f} GB，请求 {TOTAL_REQ[0]:,}", flush=True)
    print(S.round(3).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
