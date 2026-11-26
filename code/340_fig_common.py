#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
340 — 绘图公共模块（Figure 1–3, S1–S5 共用）。

职责
----
1. 全稿统一的视觉编码（Okabe–Ito）：多能 蓝、其他 深灰、LCL 橙；可读实心 / 边缘半填充 /
   不可读空心 / 未评估灰 ×；r_s_nz<=0 空心三角画在 0。
2. 版面：174 mm 全版心宽（FIG “174 mm and not higher than 234 mm”），Arial 8 pt
   （FIG “about 2-3 mm (8-12 pt)”），面板字母 9 pt 粗体，线宽 >=0.5 pt（FIG 下限 0.3 pt）。
   图按 17.4 cm 原尺寸插入，不做缩放，因此图上 pt 即最终 pt。
3. 每个图脚本开头调用 spot_check()：随机 3 个拷贝，从源 bigWig 重算 pa5，与表里同 unit_id
   的值比对。这是 stride-24 错位两次之后的硬规则。
4. 输出：FigX.pdf（矢量）+ FigX.png（300 dpi）。

单位换算：1 mm = 1/25.4 in。
"""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
FIGDIR = os.path.join(PAPER, "figures", "generated")
TMP = os.path.join(PAPER, "tmp")
MM = 1 / 25.4

BLUE = "#0072B2"
GREY = "#555555"
ORANGE = "#E69F00"
LGREY = "#BBBBBB"
RED = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"

PLURI = ["H1", "H9", "WTC-11", "GM27730"]
PLURISET = set(PLURI)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "lines.linewidth": 0.75,
    "lines.markersize": 3,
    "errorbar.capsize": 1.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "figure.dpi": 150,
})


def colour(sample):
    if sample in PLURISET:
        return BLUE
    return ORANGE if sample in LCL else GREY


# 数据表里的 sample 是 compendium 标识符，作图时按标准细胞系名显示
DISP = {"Hek293T": "HEK293T", "Panc1": "PANC-1", "Hap1": "HAP1"}


def disp(names):
    if isinstance(names, str):
        return DISP.get(names, names)
    return [DISP.get(n, n) for n in names]


LCL = {"HG002", "GM12878", "HG01123", "HG01258", "HG01358", "HG01361", "HG01891",
       "HG02071", "HG02074", "HG02132", "HG02135", "HG02257", "HG02486", "HG02559",
       "HG02572", "HG02717", "HG02886", "HG03516", "HG03804", "HG03942", "HG04160",
       "HG04187", "GM28570", "GM28572", "GM25455", "GM25456"}


def panel(ax, letter, dx=-0.14, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")


def num(x, nd=2):
    """数字转字符串，负号用 U+2212（与置信区间里的减号一致，不用 ASCII 连字符）。"""
    return f"{x:.{nd}f}".replace("-", "\u2212")


def clean(ax, left=True, bottom=True, top=False, right=False):
    ax.spines["top"].set_visible(top)
    ax.spines["right"].set_visible(right)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(direction="out")


def save(fig, name, widths_mm=174, dpi=600):
    """出 PDF（矢量）+ EPS（矢量）+ PNG（600 dpi 底）。

    闸门：PNG 宽度必须 >= 4110 px，即按 17.4 cm 原尺寸插入时仍有 600 dpi（组合图）。
    FIG 规定图宽 174 mm、高 <=234 mm；按 17.4 cm 插入不缩放，图上 8 pt 即最终 8 pt。
    """
    os.makedirs(FIGDIR, exist_ok=True)
    pdf = os.path.join(FIGDIR, f"{name}.pdf")
    png = os.path.join(FIGDIR, f"{name}.png")
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.05)
    eps = os.path.join(FIGDIR, f"{name}.eps")
    try:
        fig.savefig(eps, format="eps", bbox_inches="tight", pad_inches=0.05)
    except Exception as e:                      # EPS 不支持透明时降级
        print("     EPS 失败:", e)
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  -> {pdf}")
    print(f"  -> {png}")
    size = os.path.getsize(pdf) / 1024
    print(f"     pdf {size:.0f} KB")
    import struct
    with open(png, "rb") as fh:
        w, h = struct.unpack(">II", fh.read(24)[16:24])
    eff = w / (17.4 / 2.54)          # 按 17.4 cm 插入时的有效 dpi
    print(f"     png {w}x{h}px -> {eff:.0f} dpi @17.4 cm, canvas {w*25.4/dpi:.0f} mm")
    assert w >= 4110, (f"{name}: PNG 仅 {w}px，按 17.4 cm 插入只有 {eff:.0f} dpi "
                       f"（<600），画布 {w*25.4/dpi:.0f} mm 太窄，按 174 mm 重排")


def dump_numbers(name, d):
    os.makedirs(FIGDIR, exist_ok=True)
    p = os.path.join(FIGDIR, f"{name.lower()}_numbers.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(d, fh, ensure_ascii=False, indent=1, default=float)
    print(f"  -> {p}")
    return p


# ---------------------------------------------------------------- 抽查
def spot_check(tag, n_copy=3, seed=None, samples=("H1", "Hek293T"), force=False):
    """随机 n_copy 个 L1PA2 拷贝，从源 bigWig 重算 pa5，与 t9b 表比对。

    结果写 04_tmp/fig_spotcheck.json（同一 tag 只算一次），返回 max|Δ|。
    """
    if os.environ.get("LIPA2_SKIP_SPOTCHECK") == "1":
        print(f"[spot_check] {tag}: skipped (LIPA2_SKIP_SPOTCHECK=1)")
        return 0.0
    cache_p = os.path.join(TMP, "fig_spotcheck.json")
    os.makedirs(TMP, exist_ok=True)
    cache = {}
    if os.path.exists(cache_p):
        cache = json.load(open(cache_p, encoding="utf-8"))
    if tag in cache and not force:
        return cache[tag]["max_abs_diff"]

    import importlib.util
    spec = importlib.util.spec_from_file_location("remotebw", os.path.join(HERE, "_remotebw.py"))
    rb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rb)

    units = pd.read_csv(os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz"), sep="\t")
    units = units[units.part == "L1PA_5p"].copy()
    units["grp"] = units.unit_id.str.rsplit("_", n=1).str[-1]
    units = units[units.grp == "L1PA2"].reset_index(drop=True)
    pa = pd.read_csv(os.path.join(RES, "t9b", "t9b_pa5_L1PA2_41_fixed.tsv.gz"), sep="\t").set_index("unit_id")
    man = pd.read_csv(os.path.join(PAPER, "data", "meta", "fiberseq_samples.tsv"), sep="\t").set_index("sample")

    rng = np.random.default_rng(seed if seed is not None else 20261125)
    idx = rng.choice(len(units), size=n_copy, replace=False)
    rows, diffs = [], []
    for s in samples:
        bw = rb.RemoteBigWig(man.loc[s, "percent_accessible_all"])
        for i in idx:
            t = units.iloc[i]
            v = bw.window_mean(t.chrom, int(t.start), int(t.end))
            ref = float(pa.loc[t.unit_id, s])
            rows.append(dict(tag=tag, sample=s, unit_id=t.unit_id, bigwig=v, table=ref,
                             abs_diff=abs(v - ref)))
            diffs.append(abs(v - ref))
        bw.close()
    mx = float(np.max(diffs))
    cache[tag] = dict(max_abs_diff=mx, rows=rows)
    json.dump(cache, open(cache_p, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"  spot check [{tag}] {len(rows)} 值，max|Δ| = {mx:.2e}")
    assert mx < 1e-6, f"抽查不一致 {mx}"
    return mx


# ---------------------------------------------------------------- 数据装载
def load(name, **kw):
    return pd.read_csv(os.path.join(RES, name), sep="\t", **kw)


def rel41():
    d = load("t15/t15_reliability_41.tsv", keep_default_na=False, na_values=[""])
    return d


def sample_order_by_rfull(d):
    """按 r_full 从高到低，未评估放最底（图 2a/b 共用）。"""
    a = d[np.isfinite(d.r_full)].sort_values("r_full", ascending=False)
    b = d[~np.isfinite(d.r_full)]
    return list(a["sample"]) + list(b["sample"])


def style_sample_legend(ax, ncol=1, loc="upper right", frameon=False):
    h = [Line2D([], [], marker="o", ls="", color=BLUE, ms=3, label="pluripotent"),
         Line2D([], [], marker="o", ls="", color=ORANGE, ms=3, label="lymphoblastoid"),
         Line2D([], [], marker="o", ls="", color=GREY, ms=3, label="other")]
    ax.legend(handles=h, loc=loc, ncol=ncol, frameon=frameon, handletextpad=0.2,
              borderaxespad=0.1, labelspacing=0.2)
