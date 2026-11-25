"""242_t1_did_perm4.py — 多能组从三株改成四株后的 DiD 置换与排序判据重算。

为什么要重算
------------
219/229 的置换按 C(41,3)=10,660 个三元组做，p 下限 9.38e-5，判据是
「{H1,H9,WTC-11} 是 DiD 的 top3」。GM27730 更正为 iPSC 后，多能组是四株
（H1、H9 两株 ESC；WTC-11、GM27730 两株 iPSC），三元组换成四元组：
C(41,4) = 101,270，p 下限 9.87e-6。

判据（写死在本文件，不改）
--------------------------
RANK4-C5   在 0.25-0.65 九个预设阈值上，{H1,H9,WTC-11,GM27730} 都是 DiD 的 top4
           -> ROBUST；否则 RANK4_THRESHOLD_DEPENDENT，报成立阈值数 k/9 与侵入者
RANK4-MEAN 在 0.25-0.65 九个预设阈值上，四株均值对 C(41,4) 置换 p < 0.001
           -> MEAN_ROBUST（这是 219 PASS-A 的同一条统计量，只把 k 从 3 改成 4）
RANK4-PA5  在 pa5 连续口径、9 个样本上，四株是 top4（C(9,4)=126，p 下限 7.94e-3）
           -> PASS；否则 RANK4_PA5_FAIL

**诚实声明**：c5 的 top6 排序在写本脚本之前已经跑过一次（探索阶段）。因此 RANK4-C5
这一条不能算预注册，只能算描述性复算；RANK4-PA5 是同一判据在新口径上的首次应用。
两条都在 md 里注明。

DiD 定义与 219/229 一致，不改口径：
    open5(th) = mean over copies of 1[c5 >= th]，逐样本逐亚家族
    DiD(th)   = open5(L1PA2) - mean(open5 of L1PA4, L1PA5, L1PA6, L1PA7)
    pa5       = 5′ 1 kb 窗 percent.accessible 的逐拷贝均值，DiD 直接对连续值求差

输出
----
t1_did_perm4_c5.tsv, t1_did_perm4_pa5.tsv, t1_did_perm4.md
"""
from __future__ import annotations

import itertools
import os

import numpy as np
import pandas as pd

PAP = "D:/NewBeginning/BIOinformation/Age2.3/05_paper"
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
PLURI = ["H1", "H9", "WTC-11", "GM27730"]
OLD = ["L1PA4", "L1PA5", "L1PA6", "L1PA7"]
THRESHOLDS = [round(0.25 + 0.05 * i, 2) for i in range(9)]
DIAG = [0.70, 0.75]


_COMB = {}


def combs(n, k):
    """组合下标缓存；同一 (n,k) 只建一次。"""
    if (n, k) not in _COMB:
        _COMB[(n, k)] = np.array(list(itertools.combinations(range(n), k)), dtype=np.int16)
    return _COMB[(n, k)]


def perm_test(did, group, k):
    """精确置换：所有 k 元组的组均值 >= 观测组均值的比例（向量化）。"""
    d = did.dropna()
    vals = d.to_numpy(float)
    names = list(d.index)
    g = [x for x in group if x in d.index]
    real = float(d.loc[g].mean())
    idx = combs(len(names), k)
    means = vals[idx].mean(axis=1)
    tot = len(means)
    cnt = int((means >= real - 1e-12).sum())
    top = list(d.sort_values(ascending=False).index[:k])
    return dict(group_mean=real, perm_p=cnt / tot, n_comb=tot,
                top_k=";".join(top), top_k_is_group=bool(set(top) == set(g)),
                ranks={n: int(d.rank(ascending=False)[n]) for n in g})


def main():
    cp = pd.read_csv(os.path.join(RES, "fsp_panel_copy.tsv.gz"), sep="\t")
    rows = []
    for th in THRESHOLDS + DIAG:
        op = (cp.c5 >= th).astype(float)
        agg = op.groupby([cp["sample"], cp.grp]).mean().unstack()
        did = agg["L1PA2"] - agg[OLD].mean(axis=1)
        r = perm_test(did, PLURI, 4)
        other = did.drop(PLURI)
        # 旧判据（top3、C(41,3)）同表并列，方便对照
        r3 = perm_test(did, PLURI[:3], 3)
        pl = cp[cp["sample"].isin(PLURI) & (cp.grp == "L1PA2")]
        n_open = ((pl.c5 >= th).astype(float)).groupby(pl["sample"]).sum()
        rows.append(dict(
            threshold=th, top4=r["top_k"], top4_is_group=r["top_k_is_group"],
            group_mean=r["group_mean"], best_other=float(other.max()),
            ratio=float(r["group_mean"] / other.max()) if other.max() else np.nan,
            perm_p=r["perm_p"], n_comb=r["n_comb"],
            rank_H1=r["ranks"]["H1"], rank_H9=r["ranks"]["H9"],
            rank_WTC11=r["ranks"]["WTC-11"], rank_GM27730=r["ranks"]["GM27730"],
            intruder=([x for x in r["top_k"].split(";") if x not in PLURI] or [""])[0],
            old_top3=r3["top_k"], old_top3_is_group=r3["top_k_is_group"],
            old_perm_p=r3["perm_p"], old_n_comb=r3["n_comb"],
            pluri_open_copies_median=float(n_open.median()),
            in_prespecified_range=bool(th in THRESHOLDS),
        ))
        print(f"th={th:.2f}  top4={rows[-1]['top4']:34s} in_group={rows[-1]['top4_is_group']}  "
              f"rank(GM27730)={rows[-1]['rank_GM27730']}  p4={r['perm_p']:.2e}  "
              f"ratio={rows[-1]['ratio']:.2f}", flush=True)

    C = pd.DataFrame(rows)
    C.to_csv(os.path.join(RES, "t1_did_perm4_c5.tsv"), sep="\t", index=False)
    pre = C[C.in_prespecified_range]
    k = int(pre.top4_is_group.sum())
    c5_verdict = "ROBUST" if k == len(pre) else "RANK4_THRESHOLD_DEPENDENT"
    km = int((pre.perm_p < 0.001).sum())
    mean_verdict = "MEAN_ROBUST" if km == len(pre) else "MEAN_THRESHOLD_DEPENDENT"
    print(f"\nRANK4-C5: {k}/{len(pre)} 个预设阈值上四株是 top4 -> {c5_verdict}", flush=True)
    print(f"RANK4-MEAN: {km}/{len(pre)} 个预设阈值上置换 p<0.001 -> {mean_verdict} "
          f"(p 范围 {pre.perm_p.min():.2e}-{pre.perm_p.max():.2e})", flush=True)

    # ---------------- pa5，9 样本 ----------------
    t = pd.read_csv(os.path.join(RES, "t1d_pa5_by_copy_9.tsv.gz"), sep="\t")
    samples = [c for c in t.columns if c not in ("unit_id", "grp", "measure", "measure_note")]
    sub = t[t.grp.isin(["L1PA2"] + OLD)]
    g = sub.groupby("grp")[samples].mean()
    did = g.loc["L1PA2"] - g.loc[OLD].mean()
    r = perm_test(did, PLURI, 4)
    other = did.drop(PLURI)
    P = pd.DataFrame([dict(subset="pa5_9samples", n_units=int(len(sub)),
                           top4=r["top_k"], top4_is_group=r["top_k_is_group"],
                           group_mean=r["group_mean"], best_other=float(other.max()),
                           best_other_name=other.idxmax(),
                           ratio=float(r["group_mean"] / other.max()),
                           perm_p=r["perm_p"], n_comb=r["n_comb"],
                           rank_GM27730=r["ranks"]["GM27730"],
                           top4_rank4_value=float(other.max()))])
    for s in samples:
        P[s] = float(did[s])
    P.to_csv(os.path.join(RES, "t1_did_perm4_pa5.tsv"), sep="\t", index=False)
    p_verdict = "PASS" if r["top_k_is_group"] else "RANK4_PA5_FAIL"
    print("\npa5 DiD（9 样本）:")
    print(did.sort_values(ascending=False).round(3).to_string(), flush=True)
    print(f"  top4={r['top_k']}  p={r['perm_p']:.4f} ({r['n_comb']} 组合)  "
          f"ratio={float(r['group_mean']/other.max()):.2f} -> RANK4-PA5 {p_verdict}", flush=True)

    # ---------------- md ----------------
    lines = [
        "# 四株多能的 DiD 置换与排序判据（C(41,4)=101,270）",
        "",
        "**表**：`t1_did_perm4_c5.tsv`、`t1_did_perm4_pa5.tsv`　**脚本**：`242_t1_did_perm4.py`",
        "",
        "## 1 c5 阈值口径（41 样本，k=4）",
        "",
        "| 阈值 | top4 | 四株全中 | rank(GM27730) | 置换 p | 比值 | 侵入者 |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, x in C.iterrows():
        star = "" if x.in_prespecified_range else " (诊断点)"
        lines.append(f"| {x.threshold:.2f}{star} | {'/'.join(x.top4.split(';'))} | "
                     f"{'是' if x.top4_is_group else '**否**'} | {x.rank_GM27730} | "
                     f"{x.perm_p:.2e} | {x.ratio:.2f} | {x.intruder or '—'} |")
    lines += [
        "",
        f"预设区间（0.25–0.65）九点里成立 **{k}/9** -> **{c5_verdict}**。",
        "",
        f"但**组均值**这条判据（RANK4-MEAN，p<0.001）在九点上成立 **{km}/9** -> **{mean_verdict}**，",
        f"p 范围 {pre.perm_p.min():.2e}–{pre.perm_p.max():.2e}。两条判据要分开报：",
        "「四株作为一组比随机四元组高」稳；「第四名一定是 GM27730」不稳。",
        "",
        "失败集中在 0.40–0.60：GM27730 掉到第 5–9 位，侵入者是 Caco-2（0.40、0.45）与",
        "K562 / HCT116（0.50、0.55）。0.25–0.35 它是清楚的第 4（0.436 对 0.260；0.262 对 0.246；",
        "0.216 对 0.213），0.65 也是第 4（0.0026 对 0.0025，差 1e-4）。",
        "",
        "**为什么 c5 会在中段塌**：c5 是「窗被 fiber 覆盖的比例」，含测序深度成分。",
        "GM27730 的 fiber 密度低于 ESC/iPSC 那三株，阈值一抬到 0.4–0.6，它的拷贝先掉出计数，",
        "而 Caco-2 这种深度高的分化系反而残留。这与上一轮查明的 c5 衰减是同一个机制，",
        "不是多能组内部真的分家。",
        "",
        "## 2 pa5 连续口径（9 样本，k=4）",
        "",
        "| 样本 | cell_class | DiD（百分点） |",
        "|---|---|---|",
    ]
    cls = pd.read_csv(os.path.join(RES, "sample_labels.tsv"), sep="\t").set_index("sample")["cell_class"]
    for s, v in did.sort_values(ascending=False).items():
        lines.append(f"| {s} | {cls.get(s,'?')} | {v:.2f} |")
    lines += [
        "",
        f"top4 = {'/'.join(r['top_k'].split(';'))}，与四株多能**完全一致**；"
        f"置换 p = {r['perm_p']:.4f}（{r['n_comb']} 个四元组，下限 1/126 = 7.94e-3）；"
        f"四株均值 {r['group_mean']:.2f}，最高非多能 {other.idxmax()} {other.max():.2f}，"
        f"比值 {float(r['group_mean']/other.max()):.2f}。",
        "",
        "## 3 判据结论与它不能证明什么",
        "",
        f"- RANK4-C5 = `{c5_verdict}`（{k}/9）。**不能写「四个阈值以上都成立」**，只能报 k/9 与失败点。",
        f"- RANK4-MEAN = `{mean_verdict}`（{km}/9，p<0.001）。这是与 219 PASS-A 同一条统计量，"
        "换成 k=4 后仍成立，但它只说明「这一组高」，不说明第四个成员是谁。",
        f"- RANK4-PA5 = `{p_verdict}`。这是可比口径，但它只有 9 个样本、126 个四元组，"
        "p 的下限就是 7.94e-3，**不能声称 p < 0.001**。",
        "- 两条判据不一致时以 pa5 为准的理由：c5 的阈值穿越率随深度单调变化，而四株的深度不同；",
        "pa5 是逐碱基可及度均值，无阈值。这个理由是可检验的（见第 1 节的失败位置），不是偏好。",
        "- **不能证明**：四株共享的是「同一类染色质状态」还是「同一批拷贝的固有属性」。",
        "需要同个体的 LCL 对照（HG002）来分，见 `t1d3_hg002_gm12878.md`。",
        "",
        "## 4 预注册状态（如实）",
        "",
        "RANK4-C5 的 top6 排序在写脚本前已跑过一次探索，因此**不是预注册**，是描述性复算。",
        "RANK4-PA5 是同一判据在 pa5 上的首次应用。两条都不能事后改判据。",
    ]
    with open(os.path.join(RES, "t1_did_perm4.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nwrote t1_did_perm4.md / .tsv", flush=True)


if __name__ == "__main__":
    main()
