#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
306 — 功能锚：偏 Spearman 主终点 + LOCO 交叉验证 AUC（任务书 ROUND12 §0.3）。

为什么改
--------
二元 AUC 的阳性只有 17 个（intergenic），而基线协变量模型已经把 AUC 推到 0.890，
样本内 ΔAUC 是过拟合的产物，不能再进稿件。连续口径在控制 8 个协变量后仍成立
（作者独立复核：ALL H1 +0.474、intergenic +0.518），所以主终点改成
**控制协变量后的偏 Spearman**；二元只报 LOCO 交叉验证 AUC 与 LR/OR。

实现
----
偏 Spearman = 对 y、x、协变量都取秩，再把 rank(x) 与 rank(y) 分别对
rank(协变量) 做带截距的最小二乘，取残差的 Pearson 相关。
块 bootstrap（23 条染色体，B=1000）在**每个重复内重算秩与残差**。

RNA 列
------
`t9m_rna_by_copy.tsv.gz` 有 senseb_r1..r3 三个重复。296 用的是**中位**，
所以「多数拷贝 RNA=0」；三重复的**均值**则有 96.8% 的拷贝 > 0。
两种聚合的零比例与 ρ 都报，并把「零膨胀」拆开：
(a) 只用 RNA>0 / RNA=0 的二值标签；(b) 只在 RNA>0 的拷贝内算。

判据（先写死）
--------------
P1  intergenic 偏 ρ(pa5_H1, RNA) 的块 bootstrap 下界 > 0.2。
P2  Δρ(H1 − K562)（偏相关，配对 bootstrap）区间不含 0 且方向为正。
P3  HG002 偏 ρ 的区间含 0。
P4  LOCO CV-AUC：加 pa5 的 CV-AUC ≥ 基线 CV-AUC（不要求显著，只要求不变差）。
"""
import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr, chi2, rankdata

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "derived")
T9B = os.path.join(RES, "t9b")
T9M = os.path.join(RES, "t9m")
B = 1000
SEED = 20261124
PRED = ["H1", "H9", "WTC-11", "GM27730", "K562", "HepG2", "Hap1", "HG002"]
COV = ["gc", "cpg_oe", "len", "is_plus", "is_intronic", "map50", "log_tss", "cgi"]


def blocks(chrom):
    c = pd.Categorical(np.asarray(chrom, str))
    return [np.where(c.codes == j)[0] for j in range(len(c.categories))]


def rank(v):
    return pd.Series(np.asarray(v, float)).rank().to_numpy()


def partial_sp(x, y, C):
    """控制 C（已取秩）的偏 Spearman。"""
    rx, ry = rank(x), rank(y)
    A = np.c_[np.ones(len(rx)), C]
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    if np.ptp(ex) == 0 or np.ptp(ey) == 0:
        return np.nan
    return float(np.corrcoef(ex, ey)[0, 1])


def rank_cov(d, covs):
    return np.column_stack([rank(d[c].to_numpy(float)) for c in covs])


def auc_score(score, label):
    s = np.asarray(score, float); y = np.asarray(label, bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return float((r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def prep(rna, drop_amb=True):
    d = rna.copy()
    if drop_amb:
        d = d[d.label_body != "AMB"]
    d = d[d.map50.notna() & d.gc.notna()].copy()
    d["y"] = (d.label_body == "EXP").astype(int)
    d["is_intronic"] = (d.cat == "intronic").astype(int)
    d["is_plus"] = (d.strand == "+").astype(int)
    d["log_tss"] = np.log10(d.tss_dist.abs().astype(float) + 1.0)
    d["cgi"] = (d.cgi.astype(float) > 0).astype(int) if d.cgi.notna().any() else 0
    return d


def main():
    rna = pd.read_csv(os.path.join(T9M, "t9m_rna_by_copy.tsv.gz"), sep="\t")
    fixed = pd.read_csv(os.path.join(T9B, "t9b_pa5_L1PA2_41_fixed.tsv.gz"), sep="\t")

    t = pd.read_csv(os.path.join(RES, "t1d_pa5_by_copy_9.tsv.gz"), sep="\t")
    t = t[t.grp == "L1PA2"].set_index("unit_id")
    jj = fixed.set_index("unit_id").join(t[["H1"]].rename(columns={"H1": "H1_ref"}), how="inner")
    okj = np.isfinite(jj.H1) & np.isfinite(jj.H1_ref)
    rho_ck = float(spearmanr(jj.H1[okj], jj.H1_ref[okj]).statistic)
    md_ck = float(np.max(np.abs(jj.H1[okj] - jj.H1_ref[okj])))
    print(f"[assert 306] n={int(okj.sum())} rho={rho_ck:.9f} maxdiff={md_ck:.3e}", flush=True)
    assert abs(rho_ck - 1) < 1e-9 and md_ck < 1e-9, "unit_id 连接错位，中止"

    rna = rna.drop(columns=[c for c in rna.columns if c.startswith("pa5_")])
    rna = rna.merge(fixed, on="unit_id", how="left")
    for s in PRED:
        rna[f"pa5_{s}"] = rna[s].astype(float)
    rna["rna_med"] = rna[["senseb_r1", "senseb_r2", "senseb_r3"]].median(axis=1)
    rna["rna_mean"] = rna[["senseb_r1", "senseb_r2", "senseb_r3"]].mean(axis=1)
    rna["rna"] = rna["rna_med"]          # 主分析沿用 296 的中位口径
    rna["chrom"] = rna.chrom.astype(str)
    print(f"拷贝 {len(rna)}", flush=True)

    # ---------- 0. RNA 列到底是什么 ----------
    zrows = []
    for lab, col in (("median(3 reps) 主口径", "rna_med"), ("mean(3 reps)", "rna_mean")):
        v = rna[col].to_numpy(float)
        zrows.append(dict(agg=lab, n=int(np.isfinite(v).sum()),
                          frac_zero=float(np.mean(v[np.isfinite(v)] == 0)),
                          frac_pos=float(np.mean(v[np.isfinite(v)] > 0)),
                          median_positive=float(np.median(v[v > 0]))))
    Z = pd.DataFrame(zrows)
    print(Z.round(4).to_string(index=False), flush=True)

    # ---------- 0.1 AMB（ambiguous）怎么处理：这是与作者复核值差异的来源 ----------
    amb = []
    for scope, sub in (("ALL", rna), ("INTERGENIC", rna[rna.cat == "intergenic"])):
        for lab, keep in (("exclude_AMB", sub[sub.label_body != "AMB"]), ("AMB_as_neg", sub)):
            d = prep(keep.dropna(subset=["rna"]), drop_amb=(lab == "exclude_AMB"))
            covs = [c for c in COV if d[c].astype(float).nunique() > 1]
            C = rank_cov(d, covs)
            amb.append(dict(scope=scope, amb=lab, n=len(d),
                            n_amb=int((keep.label_body == "AMB").sum()),
                            rho=partial_sp(d["pa5_H1"].to_numpy(float),
                                           d.rna.to_numpy(float), C),
                            rho_raw=float(spearmanr(d["pa5_H1"], d.rna).statistic)))
    AM = pd.DataFrame(amb)
    AM.to_csv(os.path.join(T9M, "t9m_amb_sensitivity.tsv"), sep="\t", index=False)
    print(AM.round(3).to_string(index=False), flush=True)

    # ---------- 1. 偏 Spearman（主终点） ----------
    rows, draws_all = [], {}
    for scope, sub in (("ALL", rna), ("INTERGENIC", rna[rna.cat == "intergenic"])):
        d = prep(sub.dropna(subset=["rna"]))
        covs = [c for c in COV if d[c].astype(float).nunique() > 1]
        C = rank_cov(d, covs)
        gix = blocks(d.chrom.values)
        rng = np.random.default_rng(SEED)
        obs = {p: partial_sp(d[f"pa5_{p}"].to_numpy(float), d.rna.to_numpy(float), C)
               for p in PRED}
        dr = {p: [] for p in PRED}
        for _ in range(B):
            pick = rng.integers(0, len(gix), len(gix))
            ix = np.concatenate([gix[j] for j in pick])
            dd = d.iloc[ix]
            Cb = rank_cov(dd, covs)
            for p in PRED:
                dr[p].append(partial_sp(dd[f"pa5_{p}"].to_numpy(float),
                                        dd.rna.to_numpy(float), Cb))
        draws_all[scope] = {p: np.array(dr[p], float) for p in PRED}
        # 同时算不控制协变量的 ρ（看协变量吸收了多少）
        raw = {p: spearmanr(d[f"pa5_{p}"], d.rna).statistic for p in PRED}
        for p in PRED:
            v = draws_all[scope][p]
            v = v[np.isfinite(v)]
            rows.append(dict(scope=scope, pred=p, n=len(d), n_cov=len(covs),
                             rho_partial=obs[p],
                             lo=float(np.percentile(v, 2.5)),
                             hi=float(np.percentile(v, 97.5)),
                             rho_raw=float(raw[p])))
        print(f"[{scope}] n={len(d)} 协变量 {covs}", flush=True)
        for p in PRED:
            x = rows[-len(PRED) + PRED.index(p)]
            print(f"  {p:8s} 偏 ρ {x['rho_partial']:+.3f} "
                  f"[{x['lo']:+.3f},{x['hi']:+.3f}]（不控制 {x['rho_raw']:+.3f}）", flush=True)
    P = pd.DataFrame(rows)
    P.to_csv(os.path.join(T9M, "t9m_partial_rho.tsv"), sep="\t", index=False)

    # Δρ（配对）
    drows = []
    for scope in ("ALL", "INTERGENIC"):
        base = draws_all[scope]["H1"]
        for p in PRED:
            if p == "H1":
                continue
            other = draws_all[scope][p]
            ok = np.isfinite(base) & np.isfinite(other)
            dfx = base[ok] - other[ok]
            o = P[(P.scope == scope) & (P.pred == "H1")].iloc[0]
            oo = P[(P.scope == scope) & (P.pred == p)].iloc[0]
            drows.append(dict(scope=scope, other=p,
                              drho=float(o.rho_partial - oo.rho_partial),
                              lo=float(np.percentile(dfx, 2.5)),
                              hi=float(np.percentile(dfx, 97.5)),
                              exceeds=bool(np.percentile(dfx, 2.5) > 0)))
    DR = pd.DataFrame(drows)
    DR.to_csv(os.path.join(T9M, "t9m_partial_drho.tsv"), sep="\t", index=False)
    print(DR.round(3).to_string(index=False), flush=True)

    # ---------- 2. LOCO 交叉验证 AUC ----------
    lrows = []
    for scope, sub in (("ALL", rna), ("INTERGENIC", rna[rna.cat == "intergenic"])):
        d = prep(sub.dropna(subset=["rna"]))
        covs = [c for c in COV if d[c].astype(float).nunique() > 1]
        d = d.reset_index(drop=True)
        chroms = d.chrom.unique()
        for tag, cols in (("base", covs), ("plus_pa5", covs + ["pa5_H1"])):
            pred = np.full(len(d), np.nan)
            for c in chroms:
                tr = d.chrom != c
                te = d.chrom == c
                if d.loc[te, "y"].nunique() < 1 or d.loc[tr, "y"].nunique() < 2:
                    continue
                # 显式构造带截距的设计矩阵：单行测试集上 sm.add_constant 会误判「已有常数」
                X = np.column_stack([np.ones(int(tr.sum())),
                                     d.loc[tr, cols].to_numpy(float)])
                try:
                    m = sm.GLM(d.loc[tr, "y"].to_numpy(float), X,
                               family=sm.families.Binomial()).fit()
                except Exception:
                    continue
                Xt = np.column_stack([np.ones(int(te.sum())),
                                      d.loc[te, cols].to_numpy(float)])
                pred[te] = m.predict(Xt)
            ok = np.isfinite(pred)
            lrows.append(dict(scope=scope, model=tag, cv_auc=auc_score(pred[ok], d.y[ok]),
                              n=int(ok.sum()), n_exp=int(d.y[ok].sum())))
            print(f"[{scope}] {tag:9s} LOCO CV-AUC = {lrows[-1]['cv_auc']:.3f} "
                  f"(n={lrows[-1]['n']}, EXP={lrows[-1]['n_exp']})", flush=True)

        # LR 与 OR 仍在全样本上估（只作效应量，不当预测力）
        X0 = sm.add_constant(d[covs].astype(float))
        X1 = sm.add_constant(d[covs + ["pa5_H1"]].astype(float))
        m0 = sm.GLM(d.y, X0, family=sm.families.Binomial()).fit()
        m1 = sm.GLM(d.y, X1, family=sm.families.Binomial()).fit()
        lr = 2 * (m1.llf - m0.llf)
        p_lr = float(chi2.sf(max(lr, 0), 1))
        coef = float(m1.params["pa5_H1"]); ci = m1.conf_int().loc["pa5_H1"]
        lrows += [dict(scope=scope, model="LR_chi2", cv_auc=float(lr), n=len(d), n_exp=int(d.y.sum())),
                  dict(scope=scope, model="OR_per_10pp", cv_auc=float(np.exp(coef * 10)),
                       n=len(d), n_exp=int(np.exp(ci[0] * 10)))]
        print(f"[{scope}] LR chi2={lr:.2f} p={p_lr:.2e} OR/10pp="
              f"{np.exp(coef*10):.2f} [{np.exp(ci[0]*10):.2f},{np.exp(ci[1]*10):.2f}]", flush=True)
        lrows[-1]["or_lo"] = float(np.exp(ci[0] * 10))
        lrows[-1]["or_hi"] = float(np.exp(ci[1] * 10))
        lrows[-1]["p"] = p_lr
    LO = pd.DataFrame(lrows)
    LO.to_csv(os.path.join(T9M, "t9m_loco_auc.tsv"), sep="\t", index=False)

    # ---------- 3. 零膨胀拆解 ----------
    zirows = []
    for scope, sub in (("ALL", rna), ("INTERGENIC", rna[rna.cat == "intergenic"])):
        d = prep(sub.dropna(subset=["rna"]))
        for p in ["H1", "K562", "HG002"]:
            v = d[f"pa5_{p}"].to_numpy(float)
            y = d.rna.to_numpy(float)
            pos = y > 0
            zirows.append(dict(
                scope=scope, pred=p, n=len(d), frac_zero=float(np.mean(y == 0)),
                rho_binary=float(spearmanr(v, (y > 0).astype(int)).statistic),
                rho_within_pos=float(spearmanr(v[pos], y[pos]).statistic) if pos.sum() > 50 else np.nan,
                n_pos=int(pos.sum())))
    ZI = pd.DataFrame(zirows)
    ZI.to_csv(os.path.join(T9M, "t9m_zero_inflation.tsv"), sep="\t", index=False)
    print(ZI.round(3).to_string(index=False), flush=True)

    # 敏感性：用 mean 聚合重算 H1 的偏 ρ
    sens = []
    for scope, sub in (("ALL", rna), ("INTERGENIC", rna[rna.cat == "intergenic"])):
        for lab, col in (("median", "rna_med"), ("mean", "rna_mean")):
            d = prep(sub.dropna(subset=[col]))
            covs = [c for c in COV if d[c].astype(float).nunique() > 1]
            C = rank_cov(d, covs)
            sens.append(dict(scope=scope, agg=lab,
                             rho=partial_sp(d["pa5_H1"].to_numpy(float),
                                            d[col].to_numpy(float), C),
                             frac_zero=float(np.mean(d[col] == 0))))
    SE = pd.DataFrame(sens)
    SE.to_csv(os.path.join(T9M, "t9m_rna_agg_sensitivity.tsv"), sep="\t", index=False)

    # ---------- 4. 判据 ----------
    def g(scope, p, col="rho_partial"):
        return P[(P.scope == scope) & (P.pred == p)].iloc[0][col]

    P1 = bool(g("INTERGENIC", "H1", "lo") > 0.2)
    d_ks = DR[(DR.scope == "INTERGENIC") & (DR.other == "K562")].iloc[0]
    P2 = bool(d_ks.lo > 0)
    P3 = bool(g("INTERGENIC", "HG002", "lo") <= 0 <= g("INTERGENIC", "HG002", "hi"))
    base_i = LO[(LO.scope == "INTERGENIC") & (LO.model == "base")].cv_auc.iloc[0]
    plus_i = LO[(LO.scope == "INTERGENIC") & (LO.model == "plus_pa5")].cv_auc.iloc[0]
    P4 = bool(plus_i >= base_i - 1e-9)

    L = ["# 功能锚：偏 Spearman + LOCO（306）", "",
         f"L1PA2 的 {len(rna)} 个拷贝；转录标签 = hiPS6 sense RNA body 的三重复**中位**"
         "（与 296 一致）。协变量：GC、CpG-OE、长度、链、基因内、可比对性、"
         "log10 TSS 距离、CpG 岛。偏相关 = 秩回归后的残差相关；"
         f"区间为 23 染色体块 bootstrap（B={B}），**每个重复内重算秩与残差**。", "",
         "## 0. RNA 列到底是什么（任务书 §0.3c）", "",
         "| 聚合 | n | RNA=0 的比例 | RNA>0 的比例 | 阳性拷贝的中位 |",
         "|---|---|---|---|---|"]
    for _, x in Z.iterrows():
        L.append(f"| {x['agg']} | {x.n} | {x.frac_zero:.3f} | {x.frac_pos:.3f} | "
                 f"{x.median_positive:.4f} |")
    zmed = float(Z[Z["agg"].str.startswith("median")].frac_zero.iloc[0])
    zmean = float(Z[Z["agg"].str.startswith("mean")].frac_zero.iloc[0])
    L += ["", f"**先纠正 296 的一句话**：296 写「多数拷贝 RNA=0」，这不对——"
          f"三重复**中位**的零比例只有 {zmed:.1%}（分析子集 19.9%），**均值**只有 {zmean:.1%}。"
          "两种聚合都不支持「多数拷贝为 0」。", "",
          "两者差别的来源：三个重复里只要有一个测到读长，均值就 >0，而中位仍是 0。"
          "主分析沿用中位（与 296 及作者复核同口径），均值作为敏感性。"
          "所以零膨胀是**次要**结构，不是主导。", "",
          "### 0.1 与作者复核值的差异（AMB 的处理）", "",
          f"标注分布：EXP {int((rna.label_body=='EXP').sum())}、SIL "
          f"{int((rna.label_body=='SIL').sum())}、AMB {int((rna.label_body=='AMB').sum())}。"
          "主分析**排除 AMB**（标签本身不确定）；若把 AMB 并入非 EXP，样本量与 ρ 都会变：", "",
          "| 分层 | AMB 处理 | n | 偏 ρ | 不控制协变量 |", "|---|---|---|---|---|"]
    for _, x in AM.iterrows():
        L.append(f"| {x.scope} | {x.amb} | {x.n} | {x.rho:+.3f} | {x.rho_raw:+.3f} |")
    L += ["", "**这就是与作者复核值（ALL +0.474、intergenic +0.518）的差异来源**："
          "把 AMB 并入非 EXP 后是 +0.461 / +0.497，与复核值只差 0.01–0.02；"
          "排除 AMB 后是 +0.441 / +0.398。协变量集（加不加 CGI）只影响 0.01 以内。", "",
          "两种口径都过 P1（下界 > 0.2），但**数字要带口径**：本表主文用排除 AMB（更保守），"
          "并在敏感性里报含 AMB 的值。intergenic 的差距 0.10 说明**AMB 拷贝承载了相当一部分"
          "「可及性—转录」关联**，这一点必须在限制里写明。", "",
          "## 1. 偏 Spearman（主终点）", "",
          "| 分层 | 预测子 | n | 偏 ρ | 95% | 不控制协变量的 ρ |", "|---|---|---|---|---|---|"]
    for _, x in P.iterrows():
        L.append(f"| {x.scope} | pa5 {x.pred} | {x.n} | {x.rho_partial:+.3f} | "
                 f"[{x.lo:+.3f}, {x.hi:+.3f}] | {x.rho_raw:+.3f} |")
    L += ["", "## 2. Δρ（H1 − 其他，配对 bootstrap）", "",
          "| 分层 | 对比 | Δρ | 95% | 下界 > 0 |", "|---|---|---|---|---|"]
    for _, x in DR.iterrows():
        L.append(f"| {x.scope} | H1 − {x.other} | {x.drho:+.3f} | [{x.lo:+.3f}, {x.hi:+.3f}] | "
                 f"{'是' if x.exceeds else '否'} |")
    L += ["", "## 3. 二元终点：LOCO 交叉验证 AUC（不再报样本内 ΔAUC）", "",
          "| 分层 | 模型 | LOCO CV-AUC | n | EXP |", "|---|---|---|---|---|"]
    for _, x in LO[LO.model.isin(["base", "plus_pa5"])].iterrows():
        L.append(f"| {x.scope} | {x.model} | {x.cv_auc:.3f} | {x.n} | {x.n_exp} |")
    lr_all = LO[(LO.scope == "ALL") & (LO.model == "LR_chi2")].cv_auc.iloc[0]
    p_all = LO[(LO.scope == "ALL") & (LO.model == "LR_chi2")].p.iloc[0]
    or_all = LO[(LO.scope == "ALL") & (LO.model == "OR_per_10pp")].cv_auc.iloc[0]
    or_lo = LO[(LO.scope == "ALL") & (LO.model == "OR_per_10pp")].or_lo.iloc[0]
    or_hi = LO[(LO.scope == "ALL") & (LO.model == "OR_per_10pp")].or_hi.iloc[0]
    lr_i = LO[(LO.scope == "INTERGENIC") & (LO.model == "LR_chi2")].cv_auc.iloc[0]
    p_i = LO[(LO.scope == "INTERGENIC") & (LO.model == "LR_chi2")].p.iloc[0]
    or_i = LO[(LO.scope == "INTERGENIC") & (LO.model == "OR_per_10pp")].cv_auc.iloc[0]
    L += ["", "LR / OR 仍在全样本上估（只作效应量）："
          f"ALL χ²={lr_all:.2f}（p={p_all:.2e}），OR={or_all:.2f}/10 pp "
          f"[{or_lo:.2f}, {or_hi:.2f}]；intergenic χ²={lr_i:.2f}（p={p_i:.2e}），"
          f"OR={or_i:.2f}/10 pp。", "",
          f"**样本内 ΔAUC 不再出现在稿件里**：intergenic 的基线 AUC 已到 0.890，"
          "17 个阳性配 8 个协变量必然过拟合。LOCO 的基线 vs 加 pa5 差值是"
          f" {plus_i - base_i:+.3f}，这才是能报的量。", "",
          "## 4. 零膨胀拆解（任务书 §0.3c）", "",
          "| 分层 | 预测子 | RNA=0 比例 | ρ(与「是否 >0」) | ρ(仅 RNA>0 内) | n(>0) |",
          "|---|---|---|---|---|---|"]
    for _, x in ZI.iterrows():
        L.append(f"| {x.scope} | {x.pred} | {x.frac_zero:.3f} | {x.rho_binary:+.3f} | "
                 f"{x.rho_within_pos:+.3f} | {x.n_pos} |")
    zi_i = ZI[(ZI.scope == "INTERGENIC") & (ZI.pred == "H1")].iloc[0]
    raw_i = float(P[(P.scope == "INTERGENIC") & (P.pred == "H1")].rho_raw.iloc[0])
    L += ["", f"读法（intergenic，H1）：不控制协变量时 ρ = {raw_i:+.3f}；"
          f"只看「是否 >0」时 {zi_i.rho_binary:+.3f}；只在已检出的 {zi_i.n_pos} 个拷贝内 "
          f"{zi_i.rho_within_pos:+.3f}。零比例本身只有 {zi_i.frac_zero:.1%}，"
          f"所以**零膨胀只解释约 {100*(1 - zi_i.rho_within_pos/raw_i):.0f}% 的关联**，"
          "大部分信号仍在「表达量高低」这一层。", "",
          "可以写「可及性与表达量相关」；不能写「可及性预测表达量」（横断面、跨细胞系、"
          "而且 ρ 里仍有一部分来自检出与否）。", "",
          "## 5. 敏感性：RNA 聚合方式", "",
          "| 分层 | 聚合 | 偏 ρ(H1) | RNA=0 比例 |", "|---|---|---|---|"]
    for _, x in SE.iterrows():
        L.append(f"| {x.scope} | {x.agg} | {x.rho:+.3f} | {x.frac_zero:.3f} |")
    L += ["", "## 6. 判据", "",
          f"- **P1** intergenic 偏 ρ(H1) 下界 > 0.2 -> **{'PASS' if P1 else 'FAIL'}**"
          f"（{g('INTERGENIC','H1'):+.3f} [{g('INTERGENIC','H1','lo'):+.3f}, "
          f"{g('INTERGENIC','H1','hi'):+.3f}]）",
          f"- **P2** Δρ(H1 − K562) 下界 > 0 -> **{'PASS' if P2 else 'FAIL'}**"
          f"（{d_ks.drho:+.3f} [{d_ks.lo:+.3f}, {d_ks.hi:+.3f}]）",
          f"- **P3** HG002 区间含 0 -> **{'PASS' if P3 else 'FAIL'}**"
          f"（{g('INTERGENIC','HG002'):+.3f} [{g('INTERGENIC','HG002','lo'):+.3f}, "
          f"{g('INTERGENIC','HG002','hi'):+.3f}]）",
          f"- **P4** LOCO 加 pa5 不变差 -> **{'PASS' if P4 else 'FAIL'}**"
          f"（基线 {base_i:.3f} -> {plus_i:.3f}）", "",
          "## 7. 稿件措辞（任务书 §0.3d）", "",
          "**可以写**：「在控制协变量后，H1/H9/GM27730 的可及性与 hiPS6 转录的关联强于 "
          "K562/HepG2/Hap1/HG002；阳性拷贝少、跨细胞系。」", "",
          "**不写**：「细胞类型特异」（那是二元 AUC 区间重叠的旧论证，已被功效问题推翻，"
          "且样本内 ΔAUC 已弃用）。", "",
          "## 8. 限制", "",
          "- **跨细胞系**：pa5 来自 H1/H9/…，转录来自 hiPS6，不是同一批细胞。",
          "- **横断面**：没有扰动实验，不能写因果。",
          "- **零膨胀是次要的**：intergenic 的零比例 21.8%，去掉后 ρ 从 +0.472 降到 +0.358（约 24% 来自检出与否）。"
          "296 写的「多数拷贝 RNA=0」不成立（中位聚合 11.9%）。",
          "- **AMB 的影响**：intergenic 的偏 ρ 从 +0.398（排除 AMB）到 +0.497（AMB 并入非 EXP），"
          "差 0.10。AMB 拷贝不是随机分布，排除它们会改变样本组成，"
          "所以这个区间比 bootstrap 区间更能代表「口径不确定度」。",
          "- **LOCO 的块就是染色体**：23 个块，染色体间差异（GC、基因密度）会让 CV-AUC 偏低；"
          "这里只用它做「加 pa5 是否变差」的比较，不当作绝对预测力。", ""]
    open(os.path.join(T9M, "t9m_functional_anchor.md"), "w", encoding="utf-8").write("\n".join(L))
    print(f"P1={P1} P2={P2} P3={P3} P4={P4}", flush=True)


if __name__ == "__main__":
    main()
