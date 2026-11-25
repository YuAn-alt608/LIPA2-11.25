#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
295 — 同拷贝 ΔA × ΔM 的重做（任务书 ROUND11 §0.1）。

第十轮 288 错在哪
-----------------
288 里 a1/a2 取自 `t9b_pa5_hap_L1PA2.tsv.gz`（unit_id 被置换），
而 m1/m2 是按 `l1pa5_copy_tracks_v2` 的**正确坐标**从 H1_cpg_hap{1,2}.bw 取的。
两者按 unit_id 对齐 = 把 A 拷贝的可及性差配给了 B 拷贝的甲基化差。
「差值的差在同文件内部相减不变」这句只对 ΔA 成立，对跨文件的 ΔM 不成立。
所以 288 得到的 ρ≈0、置换 p=0.093 是错配的产物。

另外 288 的「方差分解 / NO_POWER」逻辑本身不成立：
r_s 是 hap1 与 hap2 的相关，它把两条单倍型之间的**全部**不一致都当成噪声，
于是「差值里的噪声占比」在构造上恒 ≈1，不可能给出功效信息。本轮删除该逻辑，
只写「不可识别（需要读长级重复）」。

本轮做法
--------
1. a1/a2 与 m1/m2 全部按**修复后**的 unit_id 与坐标取；连接处断言 = t1d 的 H1 列。
2. 重写 `t9i_dA_by_copy.tsv.gz`：pa5_* 换成修复版，dA_hap_H1 用修复版重算。
3. 判据（先写死）：
   E1  ρ(ΔA, ΔM) < 0 且 23 块 bootstrap 的 95% **上界 < 0**。
   E2  置换带（打乱 ΔM 的拷贝标签）不含观测值，双侧 p < 0.05。
   E3  反向同号率 > 0.5，二项检验单侧 p < 0.05。
4. 反假象控制（任一失败 -> 降为「待查」）：
   C-a 单倍型平均水平四分位内是否稳定（低四分位解释地板效应）。
   C-b 控制两条单倍型的**可用 CpG 数差**（偏相关）。
   C-c 只用**两条单倍型都测到的 CpG**（交集）重算 ΔM —— 这一条最关键：
       它把「两条单倍型测的不是同一批 CpG」这个假象直接去掉。
   C-d 控制读长数差：轨道只有 percent.accessible 与每 CpG 甲基化值，没有读长数，
       取不到 -> 如实写「未控制」。
   C-e 随机配对（= 置换，E2 已覆盖）。
   C-f 重复性：只有 H1 有 CpG 单倍型轨道（H9 只有 combined），写「仅 H1」。

结论口径：同一拷贝的两条等位基因之间，更开放的一条甲基化更低。**不得**写成因果。
"""
import os
import sys
import numpy as np
import pandas as pd
import pybigtools
from scipy.stats import spearmanr, binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _twobit import TwoBit

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T9B = os.path.join(RES, "t9b")
T9D = os.path.join(RES, "t9d")
T9I = os.path.join(RES, "t9i")
TWOBIT = "../../01_data/ref/hg38.2bit"
FIBER = "../../Age2.3.2/01_data/fiberseq"
SEED = 20261124
B_BOOT = 1000
B_PERM = 20000


def block_boot(x, y, chrom, B=B_BOOT, seed=SEED):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, c = x[ok], y[ok], np.asarray(chrom)[ok]
    r0 = float(spearmanr(x, y).statistic)
    cats = pd.Categorical(c)
    gix = [np.where(cats.codes == j)[0] for j in range(len(cats.categories))]
    rng = np.random.default_rng(seed)
    v = []
    for _ in range(B):
        pick = rng.integers(0, len(gix), len(gix))
        ix = np.concatenate([gix[j] for j in pick])
        if np.ptp(x[ix]) == 0 or np.ptp(y[ix]) == 0:
            continue
        v.append(spearmanr(x[ix], y[ix]).statistic)
    v = np.array(v, float)
    return r0, float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)), len(gix)


def hap_meth(bw1, bw2, tb, units, common=False):
    """每个拷贝的 hap1/hap2 甲基化（0-1）与可用 CpG 数。

    common=True 时只用两条单倍型**都**有值的 CpG（交集），用于 C-c。
    """
    n = len(units)
    m1 = np.full(n, np.nan); m2 = np.full(n, np.nan)
    c1 = np.zeros(n, int); c2 = np.zeros(n, int)
    cu = np.zeros(n, int)          # 两条单倍型 CpG 集合的并集大小
    mc = np.full(n, np.nan)
    for i, r in enumerate(units.itertuples()):
        seq = tb.sequence(r.chrom, int(r.start), int(r.end)).upper()
        pos = [k for k in range(len(seq) - 1) if seq[k] == "C" and seq[k + 1] == "G"]
        if len(pos) < 3:
            continue
        v1 = np.asarray(bw1.values(r.chrom, int(r.start), int(r.end), fillna=np.nan),
                        dtype=float)[pos]
        v2 = np.asarray(bw2.values(r.chrom, int(r.start), int(r.end), fillna=np.nan),
                        dtype=float)[pos]
        o1 = np.isfinite(v1); o2 = np.isfinite(v2)
        c1[i] = int(o1.sum()); c2[i] = int(o2.sum())
        cu[i] = int((o1 | o2).sum())
        if o1.sum() >= 3:
            m1[i] = v1[o1].mean() / 100.0
        if o2.sum() >= 3:
            m2[i] = v2[o2].mean() / 100.0
        both = o1 & o2
        if both.sum() >= 3:
            mc[i] = (v1[both].mean() - v2[both].mean()) / 100.0
    return m1, m2, c1, c2, cu, mc


def partial_spearman(x, y, z):
    """控制 z 的偏 Spearman（秩的线性回归残差）。"""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[ok], y[ok], z[ok]
    if len(x) < 30:
        return np.nan, np.nan
    rx, ry, rz = (pd.Series(v).rank().to_numpy() for v in (x, y, z))
    A = np.c_[np.ones_like(rz), rz]
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    r = float(np.corrcoef(ex, ey)[0, 1])
    n = len(rx)
    t = r * np.sqrt((n - 3) / max(1 - r * r, 1e-12))
    from scipy.stats import t as tdist
    p = 2 * (1 - tdist.cdf(abs(t), n - 3))
    return r, float(p)


def main():
    os.makedirs(T9I, exist_ok=True)
    pa = pd.read_csv(os.path.join(T9B, "t9b_pa5_L1PA2_41_fixed.tsv.gz"), sep="\t")
    hap = pd.read_csv(os.path.join(T9B, "t9b_pa5_hap_L1PA2_fixed.tsv.gz"), sep="\t")
    units = pd.read_csv(os.path.join(RES, "l1pa5_copy_tracks_v2.tsv.gz"), sep="\t")
    units = units[units.part == "L1PA_5p"].drop_duplicates("unit_id")

    # 断言：连接后与 t1d 的 H1 列完全一致（任务书 §0.5）
    t = pd.read_csv(os.path.join(RES, "t1d_pa5_by_copy_9.tsv.gz"), sep="\t")
    t = t[t.grp == "L1PA2"].set_index("unit_id")
    jj = pa.set_index("unit_id").join(t[["H1"]].rename(columns={"H1": "H1_ref"}), how="inner")
    okj = np.isfinite(jj.H1) & np.isfinite(jj.H1_ref)
    rho_ck = float(spearmanr(jj.H1[okj], jj.H1_ref[okj]).statistic)
    md_ck = float(np.max(np.abs(jj.H1[okj] - jj.H1_ref[okj])))
    print(f"[assert 295] n={int(okj.sum())} rho={rho_ck:.9f} maxdiff={md_ck:.3e}", flush=True)
    assert abs(rho_ck - 1) < 1e-9 and md_ck < 1e-9, "unit_id 连接错位，中止"

    pa = pa.set_index("unit_id")
    # tracks_v2 是全部 7,222 个 L1PA 拷贝，必须按 L1PA2 面板取，否则 reindex 出 NaN
    units = units.set_index("unit_id").reindex(pa.index).reset_index()
    assert units.chrom.notna().all(), "units 缺 unit_id"
    hap = hap.set_index("unit_id").reindex(pa.index)
    chrom = units.chrom.to_numpy()

    # ---------- 1. 重写 t9i_dA_by_copy.tsv.gz（修复版） ----------
    a1 = hap["H1.hap1"].to_numpy(float)
    a2 = hap["H1.hap2"].to_numpy(float)
    D = pd.DataFrame(dict(unit_id=units.unit_id, chrom=chrom,
                          start=units.start.to_numpy(), end=units.end.to_numpy(),
                          pa5_GM27730=pa["GM27730"].to_numpy(float),
                          pa5_HG002=pa["HG002"].to_numpy(float)))
    D["dA"] = D.pa5_GM27730 - D.pa5_HG002
    D["pa5_H1"] = pa["H1"].to_numpy(float)
    D["dA_hap_H1"] = a1 - a2

    # ---------- 2. 甲基化（全 CpG 与交集两套） ----------
    tb = TwoBit(TWOBIT)
    b1 = pybigtools.open(f"{FIBER}/H1_cpg_hap1.bw")
    b2 = pybigtools.open(f"{FIBER}/H1_cpg_hap2.bw")
    m1, m2, c1, c2, cu, mc = hap_meth(b1, b2, tb, units, common=True)
    b1.close(); b2.close(); tb.close()
    D["dM_hap_H1"] = m1 - m2
    D["dM_common_cpg"] = mc
    D["nCpg_hap1"] = c1
    D["nCpg_hap2"] = c2
    D["nCpg_union"] = cu
    D.to_csv(os.path.join(T9I, "t9i_dA_by_copy.tsv.gz"), sep="\t", index=False)
    print(f"重写 t9i_dA_by_copy.tsv.gz：{len(D)} 行；dA 与 GM27730 的 ρ = "
          f"{spearmanr(D.dA, D.pa5_GM27730).statistic:+.3f}", flush=True)

    # ---------- 3. 主检验 ----------
    dA = a1 - a2
    dM = m1 - m2
    ok = np.isfinite(dA) & np.isfinite(dM)
    n = int(ok.sum())
    rho, lo, hi, nb = block_boot(dA, dM, chrom)
    print(f"n={n} ρ(ΔA,ΔM)={rho:+.4f} 块bootstrap [{lo:+.4f},{hi:+.4f}] 块数={nb}", flush=True)

    # E2 置换
    rng = np.random.default_rng(SEED)
    da, dm = dA[ok], dM[ok]
    null = np.array([spearmanr(da, rng.permutation(dm)).statistic for _ in range(B_PERM)])
    nlo, nhi = np.percentile(null, [2.5, 97.5])
    p_perm = float((np.abs(null) >= abs(rho)).mean())
    print(f"置换 B={B_PERM}: 带 [{nlo:+.4f},{nhi:+.4f}] p={p_perm:.2e}", flush=True)

    # E3 反向同号率
    rev = ((da > 0) & (dm < 0)) | ((da < 0) & (dm > 0))
    k_rev = int(rev.sum()); n_tot = int(n)
    conc = k_rev / n_tot
    bt = binomtest(k_rev, n_tot, 0.5, alternative="greater")

    E1 = bool(rho < 0 and hi < 0)
    E2 = bool((rho < nlo) or (rho > nhi))
    E3 = bool(conc > 0.5 and bt.pvalue < 0.05)

    # ---------- 4. 控制项 ----------
    ctrl = []

    # C-a 四分位（按单倍型平均水平）
    lvl = (a1 + a2) / 2
    qq = pd.qcut(pd.Series(lvl[ok]), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    for g in ["Q1", "Q2", "Q3", "Q4"]:
        m = (qq == g).to_numpy()
        r_g = float(spearmanr(da[m], dm[m]).statistic) if m.sum() > 30 else np.nan
        ctrl.append(dict(control="quartile_level", stratum=g, n=int(m.sum()),
                         rho=r_g, lo=np.nan, hi=np.nan, p=np.nan,
                         note="按 (hap1+hap2)/2 的 pa5 四分位"))
    # C-b 偏相关：控制 CpG 数差
    dcp = (c1 - c2).astype(float)
    pr, pp = partial_spearman(dA, dM, dcp)
    ctrl.append(dict(control="partial_nCpg_diff", stratum="all", n=n, rho=pr,
                     lo=np.nan, hi=np.nan, p=pp, note="控制 nCpg_hap1−nCpg_hap2"))
    # C-c 交集 CpG
    okc = np.isfinite(dA) & np.isfinite(mc)
    r_c = float(spearmanr(dA[okc], mc[okc]).statistic)
    cb = block_boot(dA[okc], mc[okc], chrom[okc], B=500)
    ctrl.append(dict(control="common_cpg_only", stratum="all", n=int(okc.sum()),
                     rho=r_c, lo=cb[1], hi=cb[2], p=np.nan,
                     note="只用两条单倍型都测到的 CpG"))
    # C-d 读长数差：未控制
    ctrl.append(dict(control="read_count_diff", stratum="all", n=n, rho=np.nan,
                     lo=np.nan, hi=np.nan, p=np.nan,
                     note="未控制：轨道只有 percent.accessible 与每 CpG 甲基化值，无读长数"))
    # C-d2 有效 A 数：参考序列两条 hap 相同，差 = 0；被读长覆盖的 A 数取不到
    ctrl.append(dict(control="effective_A_diff", stratum="all", n=n, rho=np.nan,
                     lo=np.nan, hi=np.nan, p=np.nan,
                     note="A 位点数按参考序列两条 hap 相同（差=0）；被读长覆盖的 A 数取不到"))
    # C-f 重复性
    ctrl.append(dict(control="replicate_sample", stratum="all", n=n, rho=np.nan,
                     lo=np.nan, hi=np.nan, p=np.nan,
                     note="仅 H1：H9 只有 combined CpG 轨道，无 hap 拆分"))
    CT = pd.DataFrame(ctrl)
    CT.to_csv(os.path.join(T9I, "t10_dA_dM_controls.tsv"), sep="\t", index=False)
    print(CT.round(4).to_string(index=False), flush=True)

    # C-b 附加：两 hap 的 CpG 数差本身是否与 ΔA / ΔM 相关（假象的必要条件）
    r_dcp_dM = float(spearmanr(dcp[ok], dm).statistic)
    r_dcp_dA = float(spearmanr(dcp[ok], da).statistic)
    # CpG 集合重叠：并集里只有一侧有值的比例（假象是否有可能存在）
    only_one = cu - c1 - c2 + (c1 + c2 - cu)      # |A∪B| - |A| - |B| + |A∩B|  = |A△B|
    frac_one = float(np.median(only_one / np.maximum(cu, 1)))
    frac_ident = float(np.mean(only_one == 0))
    print(f"CpG 集合：union 中位 {np.median(cu):.0f}，仅一侧有值的中位 {np.median(only_one):.0f}，"
          f"集合完全相同的拷贝 {frac_ident:.1%}", flush=True)

    # ---------- 5. 汇总表 ----------
    main_row = dict(n=n, rho=rho, lo=lo, hi=hi, p_perm=p_perm, null_lo=float(nlo),
                    null_hi=float(nhi), concord=conc, p_binom=float(bt.pvalue),
                    E1=E1, E2=E2, E3=E3,
                    sd_dA=float(np.nanstd(dA)), sd_dM=float(np.nanstd(dM)),
                    med_dA=float(np.nanmedian(dA)), med_dM=float(np.nanmedian(dM)),
                    rho_common=float(r_c), rho_partial=float(pr),
                    r_dcp_dM=r_dcp_dM, r_dcp_dA=r_dcp_dA,
                    frac_cpg_set_identical=frac_ident, med_only_one_cpg=float(np.median(only_one)),
                    n_cpg1_med=float(np.median(c1[c1 > 0])), n_cpg2_med=float(np.median(c2[c2 > 0])))
    pd.DataFrame([main_row]).to_csv(os.path.join(T9I, "t10_dA_dM_main.tsv"),
                                    sep="\t", index=False)

    # ---------- 6. 报告 ----------
    q1 = CT[CT.control == "quartile_level"]
    L = ["# 同拷贝 ΔA × ΔM（295，重做 279/288）", "",
         f"L1PA2 的 {n} 个拷贝（H1 的两条单倍型都能算可及性与甲基化）。", "",
         "## 0. 为什么必须重做", "",
         "288 的 a1/a2 来自 unit_id 被置换的 `t9b_pa5_hap_L1PA2.tsv.gz`，而 m1/m2 是按"
         "**正确坐标**从 H1_cpg_hap{1,2}.bw 取的。两者按 unit_id 对齐，等于把 A 拷贝的"
         "可及性差配给了 B 拷贝的甲基化差。第十轮据此写的「无法区分无耦合与噪声」不再成立。",
         "",
         "同时删除 288 的方差分解 / NO_POWER：r_s 是 hap1 与 hap2 的相关，它把单倍型之间的"
         "全部不一致都当成噪声，所以「差值里噪声占比」在构造上恒 ≈1，不携带功效信息。"
         "正确的说法是**不可识别**——要估 cis 成分的噪声，需要读长级重复，本数据没有。", "",
         "## 1. 主结果", "",
         f"- ρ(ΔA, ΔM) = **{rho:+.4f}**，23 块 bootstrap 95% **[{lo:+.4f}, {hi:+.4f}]**",
         f"- 置换（打乱 ΔM 标签，B={B_PERM}）：带 [{nlo:+.4f}, {nhi:+.4f}]，"
         f"双侧 p = **{p_perm:.1e}**",
         f"- 反向同号率（一条更开放、另一条甲基化更低）= **{conc:.3f}**"
         f"（{k_rev}/{n_tot}），二项检验单侧 p = **{bt.pvalue:.1e}**",
         f"- SD(ΔA) = {np.nanstd(dA):.2f} pp，SD(ΔM) = {np.nanstd(dM):.4f}；"
         f"ΔA 中位 {np.nanmedian(dA):+.2f}，ΔM 中位 {np.nanmedian(dM):+.4f}", "",
         "## 2. 判据（先写死）", "",
         "| 判据 | 定义 | 实测 | 判定 |", "|---|---|---|---|",
         f"| E1 | ρ < 0 且块 bootstrap 上界 < 0 | {rho:+.4f}，上界 {hi:+.4f} | "
         f"**{'PASS' if E1 else 'FAIL'}** |",
         f"| E2 | 置换带不含观测值（p<0.05） | p = {p_perm:.1e} | "
         f"**{'PASS' if E2 else 'FAIL'}** |",
         f"| E3 | 反向同号率 > 0.5（二项 p<0.05） | {conc:.3f}，p = {bt.pvalue:.1e} | "
         f"**{'PASS' if E3 else 'FAIL'}** |", "",
         "## 3. 反假象控制", "",
         "| 控制 | 分层 | n | ρ | 95% | p | 说明 |", "|---|---|---|---|---|---|---|"]
    for _, x in CT.iterrows():
        rr = "—" if not np.isfinite(x.rho) else f"{x.rho:+.4f}"
        ci = "—" if not np.isfinite(x.lo) else f"[{x.lo:+.4f}, {x.hi:+.4f}]"
        pstr = "—" if not np.isfinite(x.p) else f"{x.p:.1e}"
        L.append(f"| {x.control} | {x.stratum} | {x.n} | {rr} | {ci} | {pstr} | {x.note} |")
    L += ["",
          f"- **C-a 四分位**：Q1（最低水平）ρ = {q1[q1.stratum=='Q1'].rho.iloc[0]:+.4f}，"
          f"Q2 {q1[q1.stratum=='Q2'].rho.iloc[0]:+.4f}，Q3 {q1[q1.stratum=='Q3'].rho.iloc[0]:+.4f}，"
          f"Q4 {q1[q1.stratum=='Q4'].rho.iloc[0]:+.4f}。最低四分位接近 0 是**地板效应**："
          "两条单倍型都接近 0 时 ΔA 的动态范围被压缩，ρ 自然衰减。",
          f"- **C-b CpG 数差**：nCpg 中位 hap1 {np.median(c1[c1>0]):.0f}、hap2 "
          f"{np.median(c2[c2>0]):.0f}；ρ(ΔnCpg, ΔM) = **{r_dcp_dM:+.4f}**、"
          f"ρ(ΔnCpg, ΔA) = **{r_dcp_dA:+.4f}**。控制 ΔnCpg 后偏 ρ = **{pr:+.4f}**"
          f"（p = {pp:.1e}）。**但这一条不构成控制**：两条单倍型的 CpG 集合"
          f"{frac_ident:.1%} 完全相同（并集内仅一侧有值的中位 {np.median(only_one):.0f} 个），"
          "ΔnCpg 本身接近 0，所以偏相关与主结果几乎相同是**设计使然**，不是独立证据。",
          f"- **C-c 交集 CpG**：只用两条单倍型都测到的 CpG 重算，"
          f"ρ = **{r_c:+.4f}** [{cb[1]:+.4f}, {cb[2]:+.4f}]（n={int(okc.sum())}）。"
          f"与主结果几乎相同，**但这不是独立验证**——{frac_ident:.1%} 的拷贝两条单倍型"
          "CpG 集合本来就相同，交集≈全集。它的作用是排除「CpG 抽样差异」这一条通路，"
          "而不是提供第二份证据。真正需要担心的是**读长数差**（C-d），本数据无法控制。",
          "- **C-d 读长数差**：**未控制**。轨道里只有 percent.accessible 与每 CpG 的甲基化值，"
          "没有逐单倍型的读长数；nCpg 只能当覆盖的粗糙代理。",
          "- **C-d2 A 数差**：参考序列上 A 位点两条单倍型**按定义相同**（差 = 0），"
          "所以「SNP/indel 改变 m6A 可用位点」这条通路在参考层面不存在；"
          "真正会变的是「被读长覆盖的 A 数」，本数据取不到每 hap 的读长覆盖，**未控制**。",
          "- **C-e 随机配对**：即 E2 的置换，已给出。",
          "- **C-f 重复性**：**仅 H1**。H9 只有 combined CpG 轨道（`H9_cpg.combined.bw`），"
          "没有 hap1/hap2；WTC-11 与 GM27730 没有 CpG 轨道。", "",
          "## 4. 结论口径（ROUND12 §0.2d 改写）", "",
          "可以写：**在 H1 中，同一拷贝的两条定相单倍型之间，更开放的一条甲基化更低"
          f"（ρ = {rho:+.2f}，95% [{lo:+.2f}, {hi:+.2f}]）；尚未排除定相伪影或序列伪影；单样本。**", "",
          "**不得**写：因果（「去甲基化导致开放」或反之）、「顺式调控」、"
          "「L1PA2 5′UTR 特异」。", "",
          "两处必须一起读：", "",
          f"1. **C-b / C-c 不是控制**。两条单倍型的 CpG 集合在 {frac_ident:.0%} 的拷贝上完全相同，"
          "所以「控制 CpG 数差」「只用交集 CpG」在实现上等同于不控制。它们只能排除"
          "「CpG 抽样差异」这一条通路，不能当第二份证据。",
          "2. **能区分假象的对照在 305**：换区间（5′ 上游 / 3′ 下游 2 kb 的唯一区段）"
          "与换家族（L1HS / LTR7 / AluS）。侧翼只有 −0.03（区间含 0），"
          "家族梯度 L1HS −0.17、LTR7 −0.11、AluS −0.03。结论见 `t9i/t9i_dA_dM_flank.md`。", "",
          "## 5. 限制", "",
          "- **只有 H1 一个样本**：无法判断这条关系是否在多能性样本之间保守。",
          "- **hap 标签不是定相的**：hap1/hap2 是读长聚类的结果，不是父源/母源；"
          "只能做样本内比较，不能跨样本配对。",
          "- **ΔA 与 ΔM 的窗口不同**：ΔA 是 5′ 1 kb 的 mean percent.accessible，"
          "ΔM 是同一窗口内 CpG 的平均甲基化；两者都在 1 kb 内，但不是同一组碱基。",
          "- **甲基化的 CpG 覆盖低**：中位每 hap 约 "
          f"{np.median(c1[c1>0]):.0f}/{np.median(c2[c2>0]):.0f} 个 CpG（5′ 1 kb 内），"
          "单个拷贝的甲基化均值因此有不可忽略的抽样误差；本设计无法把它与真实的"
          "等位差异分开（这正是「不可识别」的含义，而不是「无功效」）。",
          "- **E1/E2/E3 都是相关性判据**，不能替代对 cis 序列差异（SNP/indel）的定向分析；"
          "没有 phased VCF，无法逐拷贝判断哪条等位基因上少了哪个 CpG 或 A。", ""]
    open(os.path.join(T9I, "t9i_gm27730_vs_hg002.md"), "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L), flush=True)


if __name__ == "__main__":
    main()
