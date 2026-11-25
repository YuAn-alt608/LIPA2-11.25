#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_remotebw.py — 只读必要字节的 bigWig 读取器（HTTP Range + keep-alive）。

为什么不用 pybigtools 的远程模式
--------------------------------
pybigtools 0.3 能 open(https://...)，但每次 values() 实测 3 s；
41 样本 x 1032 窗口 需要 4 万次请求，即 35 小时。瓶颈是每请求都重新握手/重读索引。
本模块：每个线程一个 requests.Session（keep-alive），染色体树与索引只读一次，
每次窗口只取命中的 data block（通常 1 个，几十 KB）。

支持的 bigWig 子集
-------------------
header / 染色体 B+ 树 / cirTree 索引 / zlib data block，
item type 1(bedGraph)、2(varStep)、3(fixedStep)。
足以读 percent.accessible / WGBS / ChIP 这类 track。

接口
----
bw = RemoteBigWig(url)
bw.chroms()                            -> {chrom: size}
bw.window_stats(chrom, s, e)           -> (sum(value*bp), covered_bp)
bw.window_mean(chrom, s, e, fill=0.0)  -> 以窗口长度为分母（pa5 口径）
bw.values(chrom, s, e, fill=0.0)       -> np.array，未覆盖处 fill
bw.close()

自检：python _remotebw.py <URL> 与 pybigtools 逐值比较。
"""
import struct
import zlib
import numpy as np
import requests
from requests.adapters import HTTPAdapter

BW_MAGIC = 0x888FFC26
CHROM_MAGIC = 0x78CA8C91
CIRTREE_MAGIC = 0x2468ACE0


class RemoteBigWig:
    def __init__(self, url, session=None, timeout=60, retries=3):
        self.url = url
        self.timeout = timeout
        self.retries = retries
        if session is None:
            session = requests.Session()
            ad = HTTPAdapter(pool_connections=8, pool_maxsize=8, max_retries=0)
            session.mount("https://", ad)
            session.mount("http://", ad)
        self.s = session
        self.n_requests = 0
        self.n_bytes = 0
        h = self._get(0, 64)
        (magic, version, zoom_levels, self.chrom_tree_off, self.full_data_off,
         self.full_index_off, field_count, defined_field_count, autosql_off,
         total_summary_off, self.uncompress_buf_size, _reserved) = struct.unpack("<IHHQQQHHQQIQ", h)
        if magic != BW_MAGIC:
            raise ValueError(f"不是 bigWig（magic={magic:#x}）：{url}")
        self.version = version
        self.zoom_levels = zoom_levels
        self._name2id, self._id2name, self._sizes = self._read_chrom_tree()
        self._root = self._read_rtree()
        self.total_summary_off = total_summary_off
        self._leaf_cache = {}

    # ---------------- 底层 ----------------
    def _get(self, start, end):
        last = None
        for _ in range(self.retries):
            try:
                r = self.s.get(self.url, headers={"Range": f"bytes={start}-{end - 1}"},
                               timeout=self.timeout)
                if r.status_code not in (200, 206):
                    raise IOError(f"HTTP {r.status_code}")
                b = r.content
                self.n_requests += 1
                self.n_bytes += len(b)
                if len(b) < end - start:
                    raise IOError(f"short read {len(b)} < {end - start}")
                return b
            except Exception as e:
                last = e
        raise IOError(f"range {start}-{end} 失败：{last}")

    def _read_chrom_tree(self):
        head = self._get(self.chrom_tree_off, self.chrom_tree_off + 32)
        magic, block_size, key_size, val_size, item_count, _ = struct.unpack("<IIIIQQ", head)
        if magic != CHROM_MAGIC:
            raise ValueError("染色体树 magic 不对")
        span = 32 + item_count * (key_size + 8) + 8192
        buf = self._get(self.chrom_tree_off, self.chrom_tree_off + span)
        name2id, id2name, sizes = {}, {}, {}

        def walk(off):
            is_leaf, _, count = struct.unpack("<BBH", buf[off:off + 4])
            p = off + 4
            for _ in range(count):
                key = buf[p:p + key_size].split(b"\x00")[0].decode()
                p += key_size
                if is_leaf:
                    cid, csize = struct.unpack("<II", buf[p:p + 8])
                    p += 8
                    name2id[key] = cid
                    id2name[cid] = key
                    sizes[key] = csize
                else:
                    (child,) = struct.unpack("<Q", buf[p:p + 8])
                    p += 8
                    walk(child - self.chrom_tree_off)
        walk(32)
        return name2id, id2name, sizes

    def _read_rtree(self):
        head = self._get(self.full_index_off, self.full_index_off + 48)
        magic, block_size, item_count, sc, sb, ec, eb, end_off, items_per_slot, _ = \
            struct.unpack("<IIQIIIIQII", head)
        if magic != CIRTREE_MAGIC:
            raise ValueError("cirTree magic 不对")
        self._block_size = block_size
        self._n_items = item_count
        return self.full_index_off + 48      # 根节点

    def _node(self, off):
        """返回 (is_leaf, items)。叶子项 (sc,sb,ec,eb,dataOff,dataSize)；内部项 (sc,sb,ec,eb,childOff)。"""
        hdr = self._get(off, off + 4)
        is_leaf, _, count = struct.unpack("<BBH", hdr)
        n = count
        if is_leaf:
            b = self._get(off, off + 4 + 32 * n)
            items, p = [], 4
            for _ in range(n):
                items.append(struct.unpack("<IIIIQQ", b[p:p + 32]))
                p += 32
            return True, items
        b = self._get(off, off + 4 + 24 * n)
        items, p = [], 4
        for _ in range(n):
            items.append(struct.unpack("<IIIIQ", b[p:p + 24]))
            p += 24
        return False, items

    def _blocks(self, chrom_id, start, end):
        out = []
        stack = [self._root]
        while stack:
            off = stack.pop()
            is_leaf, items = self._node(off)
            for it in items:
                sc, sb, ec, eb, a, b = (it + (0, 0))[:6] if is_leaf else (it + (0, 0))[:6]
                # 相交判定（跨染色体的项保守放行）
                if ec < chrom_id or sc > chrom_id:
                    continue
                if sc == chrom_id and sb >= end:
                    continue
                if ec == chrom_id and eb <= start:
                    continue
                if is_leaf:
                    out.append((a, b))
                else:
                    stack.append(it[4])
        return out

    def _decode_block(self, off, size):
        raw = self._get(off, off + size)
        if self.uncompress_buf_size:
            raw = zlib.decompress(raw)
        chrom_id, cstart, cend, step, span, typ, _, n = struct.unpack("<IIIIIBBH", raw[:24])
        p = 24
        if typ == 1:
            arr = np.frombuffer(raw[p:p + 12 * n], dtype=np.dtype([("s", "<u4"), ("e", "<u4"), ("v", "<f4")]), count=n)
            return chrom_id, arr["s"].astype(np.int64), arr["e"].astype(np.int64), arr["v"].astype(np.float64)
        if typ == 2:
            arr = np.frombuffer(raw[p:p + 8 * n], dtype=np.dtype([("s", "<u4"), ("v", "<f4")]), count=n)
            s = arr["s"].astype(np.int64)
            return chrom_id, s, s + span, arr["v"].astype(np.float64)
        if typ == 3:
            v = np.frombuffer(raw[p:p + 4 * n], dtype="<f4").astype(np.float64)
            s = cstart + np.arange(n, dtype=np.int64) * step
            return chrom_id, s, s + span, v
        raise ValueError(f"未知 item type {typ}")

    def _intervals(self, chrom, start, end):
        cid = self._name2id.get(chrom)
        if cid is None:
            return None, None, None
        ss, ee, vv = [], [], []
        for off, size in self._blocks(cid, int(start), int(end)):
            _, s, e, v = self._decode_block(off, size)
            m = (e > start) & (s < end)
            if m.any():
                ss.append(np.maximum(s[m], start)); ee.append(np.minimum(e[m], end)); vv.append(v[m])
        if not ss:
            return None, None, None
        return np.concatenate(ss), np.concatenate(ee), np.concatenate(vv)

    def window_stats(self, chrom, start, end):
        s, e, v = self._intervals(chrom, int(start), int(end))
        if s is None:
            return 0.0, 0
        w = (e - s).astype(np.float64)
        return float((w * v).sum()), int(w.sum())

    def window_mean(self, chrom, start, end, fill=0.0):
        L = int(end) - int(start)
        if L <= 0:
            return np.nan
        s, cov = self.window_stats(chrom, start, end)
        return (s + fill * (L - cov)) / L

    def values(self, chrom, start, end, fill=0.0):
        L = int(end) - int(start)
        a = np.full(L, fill, dtype=np.float64)
        s, e, v = self._intervals(chrom, int(start), int(end))
        if s is not None:
            for x, y, z in zip(s, e, v):
                a[x - start:y - start] = z
        return a

    def total_summary(self):
        """bigWig 头部自带的总体统计（bigWigInfo 的 mean/basesCovered 就是它）。"""
        if not self.total_summary_off:
            return None
        b = self._get(self.total_summary_off, self.total_summary_off + 40)
        n, mn, mx, ssum, ssq = struct.unpack("<Qdddd", b)
        return dict(bases_covered=int(n), min=float(mn), max=float(mx),
                    sum=float(ssum), sum_sq=float(ssq),
                    mean=float(ssum / n) if n else float("nan"))

    def chroms(self):
        return dict(self._sizes)

    def close(self):
        try:
            self.s.close()
        except Exception:
            pass


if __name__ == "__main__":
    import sys, time, random
    import pybigtools
    url = sys.argv[1]
    t0 = time.time()
    bw = RemoteBigWig(url)
    print("open", round(time.time() - t0, 2), "s; chroms", len(bw.chroms()))
    pb = pybigtools.open(url)
    random.seed(1)
    ch = [c for c in bw.chroms() if c in pb.chroms()]
    bad = 0
    for c in ch[:4]:
        size = bw.chroms()[c]
        for _ in range(4):
            s0 = random.randrange(0, max(1, size - 5000))
            a = bw.values(c, s0, s0 + 1000)
            b = np.asarray(pb.values(c, s0, s0 + 1000, fillna=0.0), dtype=float)
            d = float(np.nanmax(np.abs(a - b)))
            bad += d > 1e-4
            print(f"{c}:{s0} maxdiff={d:.2e} mean={a.mean():.3f}")
    print("requests", bw.n_requests, "bytes", bw.n_bytes, "mismatch", bad)
