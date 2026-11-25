#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
381 — 补齐 15 个早期样本的 hap 覆盖轨道（闭合 Table S4）。

背景
----
331 重抓 26 个样本时同时取了 (pa5, covered_bp)，但 15 个 round-13 样本的 pa5
沿用 t9b 修序版，覆盖列一律写 NaN —— 那时代没有记 cov。于是 Table S4 里这 15 行
没有 n_cov / r_s_cov，四株「由可读变不可读」的样本（HG002、Panc1、A549、THP-1）
只能给上下界。

做法
----
对这 15 个样本 x 2 hap 重新取 (pa5, covered_bp)，但**只采用 covered_bp**；
pa5 仍以已存入 t15_pa5_hap_L1PA2_41.tsv.gz 的值为准（保已发表数字不动）。为证明
新取的覆盖与旧 pa5 逐窗对齐，每轨道要求新 pa5 与表中的同列 max|diff| < 1e-6，
否则整条轨道作废并报错。

窗口坐标从 unit_id 解析（chr10_10650857_10656884_L1PA2），顺序即表中行序，
因此检查点与已存表天然对齐。

输出
----
t15/_ckpt3/<sample>.<hap>.tsv    unit_id / pa5 / cov
t15/t15_hap_cov_L1PA2_41.tsv.gz  在原 26 列之外补齐这 15x2 列（其余列原样）
t15/t15_fetch_cov15_log.tsv      每轨道耗时 / 流量 / max|pa5-表|
"""
import os
import re
import sys
import time

import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _remotebw import RemoteBigWig

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")
CKPT = os.path.join(T15, "_ckpt3")
MAN = "../../01_data/meta/fiberseq_samples.tsv"
PA = os.path.join(T15, "t15_pa5_hap_L1PA2_41.tsv.gz")
TOL = 1e-6
NTHREAD = 24
B = 1000

# 15 个无覆盖轨道的样本（round-13 批次）
TODO15 = ["K562", "HepG2", "H1", "Hap1", "Hek293T", "Jurkat", "HG002", "Panc1",
          "THP-1", "H9", "A549", "MCF-7", "Caco-2", "WTC-11", "GM27730"]

UNIT_RE = re.compile(r"^(chr[0-9XYM]+)_(\d+)_(\d+)")


def units_from_pa(pa):
    rows = []
    for u in pa.unit_id:
        m = UNIT_RE.match(u)
        if not m:
            raise ValueError("无法解析 unit_id：" + u)
        rows.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return rows


def worker(url, rows):
    """返回 (pa5, cov)，顺序与 rows 相同。"""
    try:
        bw = RemoteBigWig(url)
    except Exception as e:
        return None, str(e), 0, 0
    out = []
    try:
        for c, s, e in rows:
            try:
                sm, cov = bw.window_stats(c, s, e)
                out.append((sm / B, cov))
            except Exception:
                out.append((np.nan, 0))
    finally:
        nreq, nbyte = bw.n_requests, bw.n_bytes
        bw.close()
    return out, None, nreq, nbyte


def fetch(url, rows, tag=""):
    chunks = [rows[i::NTHREAD] for i in range(NTHREAD)]
    res = [None] * NTHREAD
    ok = [False] * NTHREAD
    errs = []
    nb = [0, 0]
    t0 = time.time()
    with ThreadPoolExecutor(NTHREAD) as ex:
        futs = {ex.submit(worker, url, ch): i for i, ch in enumerate(chunks) if ch}
        for f in as_completed(futs):
            i = futs[f]
            o, err, nreq, nbyte = f.result()
            if err:
                errs.append(err)
            res[i] = o
            ok[i] = True
            nb[0] += nreq
            nb[1] += nbyte
    if errs or not all(ok) or any(r is None for r in res):
        return None, None, dict(tag=tag, secs=round(time.time() - t0, 1), requests=nb[0],
                                mb=round(nb[1] / 1e6, 1), median=np.nan, nonzero=np.nan,
                                error="; ".join(errs[:2])[:160] or "chunk missing")

    # 按块散回原位（第 r 块 = 第 r, r+24, r+48 … 行）
    flat = [None] * len(rows)
    for r, part in enumerate(res):
        for j, x in enumerate(part):
            k = r + j * NTHREAD
            assert k < len(rows) and flat[k] is None
            flat[k] = x
    assert all(x is not None for x in flat), "散回后仍有空洞"

    pa = np.array([x[0] for x in flat], dtype=float)
    cov = np.array([x[1] for x in flat], dtype=float)
    st = dict(tag=tag, secs=round(time.time() - t0, 1), requests=nb[0],
              mb=round(nb[1] / 1e6, 1), median=round(float(np.nanmedian(pa)), 3),
              nonzero=round(float(np.nanmean(pa > 0) * 100), 1), error="")
    print(f"  {tag:26s} {len(rows):5d}w {st['secs']:6.1f}s  med {st['median']:7.3f}  "
          f"nz {st['nonzero']:5.1f}%  {st['mb']:6.1f}MB", flush=True)
    return pa, cov, st


def main():
    os.makedirs(CKPT, exist_ok=True)
    pa_ref = pd.read_csv(PA, sep="\t")
    rows = units_from_pa(pa_ref)
    man = pd.read_csv(MAN, sep="\t").set_index("sample")
    assert len(rows) == 1032, len(rows)

    log = []
    for s in TODO15:
        assert s in man.index, s
        for hap in ("hap1", "hap2"):
            f = os.path.join(CKPT, f"{s}.{hap}.tsv")
            if os.path.exists(f):
                print(f"  {s}.{hap} 已存在，跳过", flush=True)
                continue
            url = man.loc[s, f"percent_accessible_{hap}"]
            pa, cov, st = fetch(url, rows, tag=f"{s}.{hap}")
            if pa is None:
                print(f"  !! {s}.{hap} 失败：{st['error']}", flush=True)
                log.append(dict(st, sample=s, track=hap, max_pa5_diff=np.nan, pass_align=False))
                continue
            ref = pa_ref[f"{s}.{hap}"].values.astype(float)
            dmax = float(np.nanmax(np.abs(pa - ref)))
            ok = dmax < TOL
            st.update(sample=s, track=hap, max_pa5_diff=dmax, pass_align=bool(ok))
            log.append(st)
            print(f"    对齐 max|pa5-表| = {dmax:.2e}  cov>0 {np.mean(cov > 0) * 100:.1f}%",
                  flush=True)
            if not ok:
                print(f"  !! {s}.{hap} 行序不符（{dmax:.3e}），该轨道不采用", flush=True)
                continue
            d = pd.DataFrame({"unit_id": pa_ref.unit_id.values, "pa5": ref, "cov": cov})
            assert len(d) == 1032 and d.unit_id.is_unique
            d.to_csv(f, sep="\t", index=False)
        if log:
            pd.DataFrame(log).to_csv(os.path.join(T15, "t15_fetch_cov15_log.tsv"),
                                     sep="\t", index=False)

    lg = pd.DataFrame(log)
    if len(lg):
        pd.DataFrame(log).to_csv(os.path.join(T15, "t15_fetch_cov15_log.tsv"),
                                 sep="\t", index=False)
        bad = lg[~lg.pass_align.fillna(False)]
        print(f"\n完成 {int(lg.pass_align.fillna(False).sum())}/{len(lg)} 条轨道通过对齐校验")
        if len(bad):
            print("未通过：\n", bad[["sample", "track", "max_pa5_diff", "error"]].to_string())
            return 1

    # ---- 合并进 41 列覆盖表：只填这 15x2 列的 cov，其余不动 ----
    cvp = os.path.join(T15, "t15_hap_cov_L1PA2_41.tsv.gz")
    cv = pd.read_csv(cvp, sep="\t")
    assert list(cv.unit_id) == list(pa_ref.unit_id)
    filled = 0
    for s in TODO15:
        for hap in ("hap1", "hap2"):
            f = os.path.join(CKPT, f"{s}.{hap}.tsv")
            if not os.path.exists(f):
                print(f"  跳过（无检查点）：{s}.{hap}")
                continue
            d = pd.read_csv(f, sep="\t")
            assert list(d.unit_id) == list(cv.unit_id), f"{s}.{hap} 行序不符"
            cv[f"{s}.{hap}"] = d["cov"].values
            filled += 1
    cv.to_csv(cvp, sep="\t", index=False, compression="gzip")
    print(f"合并 {filled} 列覆盖轨道 -> {cvp}")
    print("仍为空的样本：",
          sorted({c.split(".")[0] for c in cv.columns if c != "unit_id"
                  and cv[c].isna().all()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
