#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
390 — r28 的 DiD：覆盖碱基分母 + 样本层不确定性。

三件事
------
1. 分母修正：L1PA2 的 pa5 改为 Σv / covered_bp（387 的池化覆盖）。参照集
   L1PA4–7 没有逐拷贝覆盖，因此覆盖分母形式只对 L1PA2 计算；参照集保持原
   口径，两版之差就是分母口径的影响。
2. 参照集定义：4,048 个 L1PA4–7 拷贝合并后的中位数（论文定义），不是四个
   亚家族中位数的平均。后者是 r27 的旧口径（+12.14），不能与正文 12.27 混用。
3. 推断：删掉拷贝层 bootstrap（那是伪重复），改为
   (a) 分层样本层 bootstrap：多能组内重采样 4 个、非多能组内重采样 7 个，
       标签固定，B=20,000，三个种子；
   (b) 整体重采样 11 个样本但保留标签、丢弃任一组为空的抽样，作对照；
   (c) 精确置换 p（C(11,4)=330）。

输出
----
t15/t15_did_r28.tsv      每样本 DiD（两口径 × 两覆盖处理）
t15/t15_did_r28.md       汇总
"""
import itertools
import os
import numpy as np
import pandas as pd

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T15 = os.path.join(RES, "t15")
T1D = os.path.join(RES, "t1d_pa5_by_copy_11.tsv.gz")
PCOV = os.path.join(T15, "t15_pool_cov_L1PA2_41.tsv.gz")
OUT = os.path.join(T15, "t15_did_r28.tsv")
MD = os.path.join(T15, "t15_did_r28.md")
B = 1000
PLURI = ["H1", "H9", "WTC-11", "GM27730"]
OTHER = ["Hek293T", "K562", "Hap1", "HepG2", "MCF-7", "HG002", "GM12878"]
OLD = ["L1PA4", "L1PA5", "L1PA6", "L1PA7"]
S = PLURI + OTHER
NBOOT = 20000
SEEDS = (0, 1, 2)


def perm_p(did, k=4):
    v = did[S].values
    real = v[:k].mean()
    comb = np.array(list(itertools.combinations(range(len(S)), k)))
    means = v[comb].mean(axis=1)
    return float((means >= real - 1e-12).sum() / len(means)), len(means)


def boot_ci_strat(did, n=NBOOT, seeds=SEEDS):
    """组内重采样、标签固定。返回每个种子的区间与三者的最大偏离。"""
    a = did[PLURI].values
    b = did[OTHER].values
    out = []
    for s in seeds:
        rng = np.random.default_rng(s)
        v = np.empty(n)
        for i in range(n):
            v[i] = (a[rng.integers(0, len(a), len(a))].mean()
                    - b[rng.integers(0, len(b), len(b))].mean())
        out.append(np.percentile(v, [2.5, 97.5]))
    out = np.array(out)
    return out, out.min(axis=0), out.max(axis=0)


def boot_ci_whole(did, n=NBOOT, seed=0):
    """整体重采样 11 个样本，保留样本身份（标签跟着样本走），丢弃任一组为空的抽样。"""
    rng = np.random.default_rng(seed)
    grp = np.array(["pluripotent"] * 4 + ["other"] * 7)
    v = did[S].values
    out = []
    for _ in range(n):
        b = rng.integers(0, len(S), len(S))
        g = grp[b]
        if (g == "pluripotent").sum() == 0 or (g == "other").sum() == 0:
            continue
        out.append(v[b][g == "pluripotent"].mean() - v[b][g == "other"].mean())
    return np.percentile(out, [2.5, 97.5]), len(out)


def main():
    t = pd.read_csv(T1D, sep="\t")
    pc = pd.read_csv(PCOV, sep="\t").set_index("unit_id")
    t = t.set_index("unit_id")
    cov = pc.reindex(t.index)
    # 分母修正（仅 L1PA2；参照集没有逐拷贝覆盖，保持原口径）
    obs = {}
    for s in S:
        p = t.loc[t.grp == "L1PA2", s].astype(float)
        c = cov.loc[t.grp == "L1PA2", s].astype(float)
        obs[s] = np.where(c > 0, p * B / np.maximum(c, 1), np.nan)
    obs = pd.DataFrame(obs, index=t.index[t.grp == "L1PA2"])
    isX = obs.index.str.startswith("chrX")
    nref = int(t.grp.isin(OLD).sum())

    def ref(how="median"):
        """4,048 个 L1PA4–7 拷贝合并后的中位数/均值（论文定义）。"""
        f = t[t.grp.isin(OLD)][S].astype(float)
        return f.median() if how == "median" else f.mean()

    def ref_by_subfamily(how="median"):
        """旧口径：四个亚家族中位数/均值的平均（仅作对照）。"""
        f = t[t.grp.isin(OLD)]
        return (f.groupby("grp")[S].median().mean() if how == "median"
                else f.groupby("grp")[S].mean().mean())

    def did_series(frame, how="median", dropX=False, refset=None):
        l1 = frame[~isX] if dropX else frame
        g1 = (l1.median() if how == "median" else l1.mean())
        g2 = ref(how) if refset is None else refset(how)
        return g1 - g2

    raw_l1 = t.loc[t.grp == "L1PA2", S].astype(float)
    variants = {
        "raw_median": did_series(raw_l1, "median"),
        "raw_mean": did_series(raw_l1, "mean"),
        "cov_median": did_series(obs, "median"),
        "cov_mean": did_series(obs, "mean"),
        "cov_median_noX": did_series(obs, "median", dropX=True),
        "raw_median_subfam": did_series(raw_l1, "median", refset=ref_by_subfamily),
        "cov_median_subfam": did_series(obs, "median", refset=ref_by_subfamily),
    }
    print(f"参照集：{nref} 个 L1PA4–7 拷贝，合并中位数（论文定义）")
    print("各口径的 DiD（多能组均值 - 非多能组均值）：")
    for k, v in variants.items():
        print(f"  {k:18s} {v[PLURI].mean() - v[OTHER].mean():+.4f}  "
              f"(pluri {v[PLURI].mean():+.2f}, other {v[OTHER].mean():+.2f})")

    rows = []
    for s in S:
        rows.append(dict(sample=s, **{k: round(float(v[s]), 3) for k, v in variants.items()}))
    d = pd.DataFrame(rows)
    d["grp"] = ["pluripotent"] * 4 + ["other"] * 7
    d.to_csv(OUT, sep="\t", index=False)

    main = variants["cov_median"]
    ci_all, ci_lo, ci_hi = boot_ci_strat(main)
    ciw, nkeep = boot_ci_whole(main)
    p, ncomb = perm_p(main)
    fixed = variants["raw_median"]
    lines = ["# r28 DiD（覆盖碱基分母，合并参照中位数）", "",
             f"参照集：{nref} 个 L1PA4–7 拷贝合并后的中位数。", "",
             "| 口径 | DiD（百分点） | 多能组 | 非多能组 |", "|---|---|---|---|"]
    for k, v in variants.items():
        lines.append(f"| {k} | {v[PLURI].mean()-v[OTHER].mean():+.3f} | "
                     f"{v[PLURI].mean():+.2f} | {v[OTHER].mean():+.2f} |")
    lines += ["",
              f"主口径 cov_median：分层样本层 bootstrap 95% CI "
              f"[{ci_lo[0]:+.2f}, {ci_hi[1]:+.2f}]（组内重采样，B={NBOOT}，"
              f"三个种子 {ci_all[0,0]:+.2f}…{ci_all[0,1]:+.2f} / "
              f"{ci_all[1,0]:+.2f}…{ci_all[1,1]:+.2f} / {ci_all[2,0]:+.2f}…{ci_all[2,1]:+.2f}）",
              f"对照：整体重采样 11 个样本、保留标签（丢弃任一组为空的抽样，保留 "
              f"{nkeep}/{NBOOT} 次）[{ciw[0]:+.2f}, {ciw[1]:+.2f}]",
              f"精确置换 p = {p:.4f}（{ncomb} 个四元组，四个多能样本是否占据前四）",
              "",
              f"分母口径：固定分母 {variants['raw_median'][PLURI].mean()-variants['raw_median'][OTHER].mean():+.3f}，"
              f"覆盖分母 {variants['cov_median'][PLURI].mean()-variants['cov_median'][OTHER].mean():+.3f}，"
              f"差 {variants['cov_median'][PLURI].mean()-variants['cov_median'][OTHER].mean()-(variants['raw_median'][PLURI].mean()-variants['raw_median'][OTHER].mean()):+.4f}"
              " 个百分点。覆盖分母只对 L1PA2 计算。",
              f"旧参照口径（四亚家族中位数平均，仅供对照）：{variants['raw_median_subfam'][PLURI].mean()-variants['raw_median_subfam'][OTHER].mean():+.3f}；"
              f"覆盖分母 {variants['cov_median_subfam'][PLURI].mean()-variants['cov_median_subfam'][OTHER].mean():+.3f}。",
              "",
              "| 样本 | raw_median | cov_median | cov_median_noX |", "|---|---|---|---|"]
    for _, r in d.iterrows():
        lines.append(f"| {r['sample']} | {r.raw_median:+.2f} | {r.cov_median:+.2f} | "
                     f"{r.cov_median_noX:+.2f} |")
    open(MD, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\n写出", OUT, MD)


if __name__ == "__main__":
    main()
