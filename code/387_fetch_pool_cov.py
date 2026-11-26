#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
387 — 重取 41 个样本的 pooled 覆盖碱基数（r28 的分母修正用）。

为什么需要
----------
pa5 = 窗内 percent.accessible 之和 / 1,000。窗内没有 FIRE 调用的碱基（包括整个
窗口无 fibre 覆盖的碱基）按 0 计入分子，分母却仍是 1,000，所以「窗口只覆盖了
300 bp」的拷贝会被压到真实值的 30%。要改用「实际覆盖的碱基数」作分母，就需要
每个样本在 5′ 1 kb 窗内的 covered_bp；t9b 的 cov 列是坏的（全 1,000），因此重取。

只读 all.percent.accessible.bw（pooled 轨道），hap 轨道的覆盖已在 t15 里。
每轨道重算的 pa5 必须与 t9b 的池化 pa5 逐窗一致（1e-9），否则该轨道作废——
这同时校验了窗口坐标与行序。

输出
----
t15/t15_pool_cov_L1PA2_41.tsv.gz   unit_id + 41 列 covered_bp（0–1000）
t15/t15_fetch_pool_cov_log.tsv     每样本耗时 / 中位覆盖 / 与 t9b pa5 的最大差
"""
import os
import sys
import time

import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _remotebw import RemoteBigWig

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")
MAN = "../../01_data/meta/fiberseq_samples.tsv"
UNITS = os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz")
POOL = os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz")
OUTF = os.path.join(T15, "t15_pool_cov_L1PA2_41.tsv.gz")
LOG = os.path.join(T15, "t15_fetch_pool_cov_log.tsv")
NTHREAD = 24
B = 1000


def load_units():
    u = pd.read_csv(UNITS, sep="\t")
    u = u[(u.part == "L1PA_5p") & (u.unit_id.str.endswith("_L1PA2"))].copy()
    u = u.drop_duplicates("unit_id")
    assert (u.end - u.start == 1000).all()
    ref = pd.read_csv(POOL, sep="\t", usecols=["unit_id"])
    m = u.set_index("unit_id").loc[ref.unit_id]
    assert list(m.index) == list(ref.unit_id)
    assert len(m) == 1032
    return m.reset_index()[["unit_id", "chrom", "start", "end"]]


def worker(url, rows):
    try:
        bw = RemoteBigWig(url)
    except Exception as e:
        return None, str(e)
    out = []
    try:
        for c, s, e in rows:
            try:
                sm, cov = bw.window_stats(c, int(s), int(e))
                out.append((sm / B, cov))
            except Exception:
                out.append((np.nan, 0))
    finally:
        bw.close()
    return out, None


def fetch(url, rows, tag=""):
    chunks = [rows[i::NTHREAD] for i in range(NTHREAD)]
    res = [None] * NTHREAD
    errs = []
    with ThreadPoolExecutor(NTHREAD) as ex:
        futs = {ex.submit(worker, url, ch): i for i, ch in enumerate(chunks) if ch}
        for f in as_completed(futs):
            i = futs[f]
            o, err = f.result()
            if err:
                errs.append(err)
            res[i] = o
    if errs or any(r is None for r in res):
        return None, " | ".join(errs[:2])[:140] or "chunk missing"
    flat = [None] * len(rows)
    for r, part in enumerate(res):
        for j, x in enumerate(part):
            flat[r + j * NTHREAD] = x
    assert all(x is not None for x in flat), "散回后仍有空洞"
    return flat, None


def main():
    units = load_units()
    rows = list(zip(units.chrom, units.start, units.end))
    man = pd.read_csv(MAN, sep="\t").set_index("sample")
    ref = pd.read_csv(POOL, sep="\t").set_index("unit_id")
    samples = [c for c in ref.columns if c != "unit_id"]
    print(f"{len(rows)} 窗口, {len(samples)} 样本", flush=True)

    cov_out = {"unit_id": units.unit_id.values}
    log = []
    t_all = time.time()
    for s in samples:
        t0 = time.time()
        flat, err = fetch(man.loc[s, "percent_accessible_all"], rows, tag=s)
        if flat is None:
            print(f"  !! {s}: {err}", flush=True)
            log.append(dict(sample=s, ok=False, secs=round(time.time() - t0, 1),
                            med_cov=np.nan, max_pa5_diff=np.nan, error=err))
            continue
        pa = np.array([x[0] for x in flat], dtype=float)
        cov = np.array([x[1] for x in flat], dtype=float)
        dmax = float(np.nanmax(np.abs(pa - ref[s].values.astype(float))))
        ok = dmax < 1e-9
        cov_out[s] = cov
        st = dict(sample=s, ok=bool(ok), secs=round(time.time() - t0, 1),
                  med_cov=float(np.median(cov)), max_pa5_diff=dmax, error="")
        log.append(st)
        print(f"  {s:9s} {st['secs']:6.1f}s  med_cov {st['med_cov']:6.1f}  "
              f"max|pa5-pool| {dmax:.2e}  {'OK' if ok else 'FAIL'}", flush=True)
        pd.DataFrame(log).to_csv(LOG, sep="\t", index=False)

    lg = pd.DataFrame(log)
    good = lg[lg.ok.fillna(False)]
    print(f"\n{len(good)}/{len(lg)} 条轨道通过；总 {round(time.time() - t_all, 1)}s", flush=True)
    if len(good) != len(lg):
        print(lg[~lg.ok.fillna(False)].to_string(index=False))
        return 1
    d = pd.DataFrame(cov_out)[["unit_id"] + samples]
    d.to_csv(OUTF, sep="\t", index=False, compression="gzip")
    print("写出", OUTF, d.shape, flush=True)
    print("总体覆盖率中位 %.3f" % (d[samples].values.ravel().mean() / 1000), flush=True)


if __name__ == "__main__":
    main()
