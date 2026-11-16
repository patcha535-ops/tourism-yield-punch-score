#!/usr/bin/env python3
"""
te_analysis.py -- every number, table and figure of "Tourism Yield and Spatial Concentration in Thailand: A Provincial
Assessment Using the Punch Score" (77 provinces, 2024), computed from provinces_2024.csv.

    python te_analysis.py provinces_2024.csv --outdir results_te [--boot 10000 --perm 9999 --seed 20260101]

Writes results.json, LaTeX table bodies (tex/), figures (fig/) and a manifest stamped with a run identifier.
"""
import argparse
import hashlib
import json
import os
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats

REGIONS = ["Bangkok", "Central", "East", "North", "Northeast", "South"]
BKK = "Bangkok"


# ----------------------------------------------------------------------------- measures (Eq. 1-9)
def punch_score(R, V):
    """Eq. 1: PS_i = 100 * (R_i / R) / (V_i / V) = 100 * y_i / ybar."""
    return 100.0 * (R / R.sum()) / (V / V.sum())


def gini(x):
    """Eq. 3, unweighted."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    return 2.0 * np.sum(np.arange(1, n + 1) * x) / (n * x.sum()) - (n + 1.0) / n


def wgini(x, w):
    """Population-weighted Gini of per-capita variable x with weights w (trapezoid rule on the weighted Lorenz curve)."""
    o = np.argsort(x, kind="stable")
    x, w = np.asarray(x, float)[o], np.asarray(w, float)[o]
    p = w / w.sum()
    inc = x * w / (x * w).sum()
    L = np.concatenate([[0.0], np.cumsum(inc)])
    P = np.concatenate([[0.0], np.cumsum(p)])
    return 1.0 - np.sum((L[1:] + L[:-1]) * (P[1:] - P[:-1]))


def theil(x):
    """Eq. 4."""
    x = np.asarray(x, float)
    s = x / x.sum()
    s = s[s > 0]
    return float(np.sum(s * np.log(len(x) * s)))


def theil_decomp(x, g):
    """Eq. 5-7: returns (T, T_between, T_within) with g an array of group labels."""
    x, g = np.asarray(x, float), np.asarray(g)
    n, X = len(x), x.sum()
    tb = tw = 0.0
    for lab in np.unique(g):
        m = g == lab
        S = x[m].sum() / X
        tb += S * np.log(S / (m.sum() / n))
        tw += S * theil(x[m])
    return tb + tw, tb, tw


def topk(x, k):
    return float(np.sort(x)[-k:].sum() / x.sum())


def shape_stats(x, g):
    T, tb, tw = theil_decomp(x, g)
    return gini(x), T, 100.0 * tb / T


# ----------------------------------------------------------------------------- spatial
def haversine(lat, lon):
    lat, lon = np.radians(lat), np.radians(lon)
    dlat, dlon = lat[:, None] - lat[None, :], lon[:, None] - lon[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(a))


def rowstd(W):
    s = W.sum(1, keepdims=True)
    return np.divide(W, s, out=np.zeros_like(W), where=s > 0)


def weights(D, kind, band=None, k=5):
    n = len(D)
    D = D.copy()
    np.fill_diagonal(D, np.inf)
    if kind == "idw":
        W = 1.0 / D
    elif kind == "knn":
        W = np.zeros((n, n))
        idx = np.argsort(D, axis=1)[:, :k]
        np.put_along_axis(W, idx, 1.0, axis=1)
    else:
        W = (D <= band).astype(float)
    return rowstd(W)


def moran(x, W):
    z = x - x.mean()
    return float(z @ W @ z / (z @ z))                 # row-standardised: W0 = n


def moran_perm(x, W, B, rng):
    I = moran(x, W)
    z = x - x.mean()
    perm = np.array([moran(z[rng.permutation(len(z))], W) for _ in range(B)])
    zsc = (I - perm.mean()) / perm.std(ddof=1)
    p = (1 + np.sum(np.abs(perm - perm.mean()) >= abs(I - perm.mean()))) / (B + 1)
    return I, float(zsc), float(p)


def local_moran(x, W, B, rng, k=5):
    """Local Moran's I_i = z_i * (W z)_i with z standardised (sum z^2 / n = 1) and conditional permutation inference: the focal
    value is held fixed and the k neighbour values are drawn without replacement from the other n-1 values (row-standardised
    weights, so the spatial lag is a mean of k values).  Pseudo-p = (1 + #{|I_i^sim| >= |I_i|}) / (B + 1).
    Quadrants: 1 High-High, 2 Low-High, 3 Low-Low, 4 High-Low (own value vs spatial lag)."""
    n = len(x)
    z = (x - x.mean()) / x.std()
    lag = W @ z
    Ii = z * lag
    quad = np.where(z >= 0, np.where(lag >= 0, 1, 4), np.where(lag >= 0, 2, 3))
    p = np.zeros(n)
    for i in range(n):
        others = np.delete(z, i)
        idx = np.argpartition(rng.random((B, n - 1)), k, axis=1)[:, :k]
        sim = z[i] * others[idx].mean(1)
        p[i] = (1 + np.sum(np.abs(sim) >= abs(Ii[i]))) / (B + 1)
    return quad, Ii, p


# ----------------------------------------------------------------------------- tables (LaTeX bodies)
def fmt(x, d=0):
    return f"{x:,.{d}f}"


def tex_table(name, outdir, spec, header, rows, long=False):
    body = "\n".join(" & ".join(str(c) for c in r) + " \\\\" for r in rows)
    if long:
        s = (f"\\begin{{longtable}}{{{spec}}}\n\\toprule\n{header} \\\\\n\\midrule\n\\endfirsthead\n\\toprule\n{header} \\\\\n\\midrule\n\\endhead\n"
             f"\\bottomrule\n\\endfoot\n{body}\n\\end{{longtable}}\n")
    else:
        s = f"\\begin{{tabular}}{{{spec}}}\n\\toprule\n{header} \\\\\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}\n"
    open(f"{outdir}/tex/{name}.tex", "w", encoding="utf-8").write(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("--outdir", default="results_te")
    ap.add_argument("--boot", type=int, default=10000)
    ap.add_argument("--perm", type=int, default=9999)
    ap.add_argument("--seed", type=int, default=20260101)
    a = ap.parse_args()
    os.makedirs(f"{a.outdir}/tex", exist_ok=True)
    os.makedirs(f"{a.outdir}/fig", exist_ok=True)
    t0 = time.time()
    df = pd.read_csv(a.data)
    assert len(df) == 77
    rng = np.random.default_rng(np.random.SeedSequence(a.seed))
    res = {}
    R, V = df.receipts.values, df.visits.values
    df["yield_thb"] = R * 1e6 / V
    df["ps"] = punch_score(R, V)
    df["f"] = df.visits_foreign / df.visits
    ybar = R.sum() * 1e6 / V.sum()
    res["national"] = dict(receipts=float(R.sum()), visits=float(V.sum()), ybar=float(ybar),
                           foreign_visit_share=float(df.visits_foreign.sum() / V.sum()),
                           foreign_receipt_share=float(df.receipts_foreign.sum() / R.sum()),
                           ybar_thai=float(df.receipts_thai.sum() * 1e6 / df.visits_thai.sum()),
                           ybar_foreign=float(df.receipts_foreign.sum() * 1e6 / df.visits_foreign.sum()))

    # ---- Table 1: regions
    rows, reg = [], {}
    for r in REGIONS + ["Thailand"]:
        d = df if r == "Thailand" else df[df.region == r]
        rr, vv = d.receipts.sum(), d.visits.sum()
        ps = 100 * (rr / R.sum()) / (vv / V.sum())
        reg[r] = dict(n=len(d), receipts=float(rr), rshare=100 * rr / R.sum(), visits=float(vv), vshare=100 * vv / V.sum(),
                      yield_thb=rr * 1e6 / vv, ps=ps, fv=100 * d.visits_foreign.sum() / vv, fr=100 * d.receipts_foreign.sum() / rr)
        x = reg[r]
        rows.append([r, x["n"], fmt(rr, 2), f"{x['rshare']:.1f}", fmt(vv), f"{x['vshare']:.1f}", fmt(x["yield_thb"]), f"{ps:.1f}",
                     f"{x['fv']:.1f}", f"{x['fr']:.1f}"])
    res["regions"] = reg
    tex_table("table1", a.outdir, "lrrrrrrrrr", r"Region & Prov. & Receipts (M THB) & Share (\%) & Visits & Share (\%) & Receipts/visit (THB) & Punch Score & Foreign visits (\%) & Foreign receipts (\%)", rows)

    # ---- provincial ranking (Table 2, appendix Table 8)
    o = df.sort_values("ps", ascending=False).reset_index(drop=True)
    o["rank"] = np.arange(1, 78)
    o["rshare"] = 100 * o.receipts / R.sum()
    res["ranking"] = o[["province", "region", "yield_thb", "ps", "f"]].to_dict("list")
    pick = list(range(10)) + list(range(72, 77))
    tex_table("table2", a.outdir, "rllrrrrr", r"Rank & Province & Region & Receipts/visit (THB) & Punch Score & Foreign visits (\%) & Occ. (\%) & Receipts share (\%)",
              [[o["rank"][i], o.province[i], o.region[i], fmt(o.yield_thb[i]), f"{o.ps[i]:.1f}", f"{100 * o.f[i]:.1f}", f"{o.occupancy[i]:.1f}", f"{o.rshare[i]:.2f}"] for i in pick])
    tex_table("table8", a.outdir, "rllrrrrr", r"Rank & Province & Region & Visits & Receipts (M THB) & Receipts/visit (THB) & Punch Score (index) & Foreign visits (\%)",
              [[o["rank"][i], o.province[i], o.region[i], fmt(o.visits[i]), fmt(o.receipts[i], 2), fmt(o.yield_thb[i]), f"{o.ps[i]:.1f}", f"{100 * o.f[i]:.1f}"] for i in range(77)], long=True)
    above = df[df.ps > 100]
    res["above100"] = dict(n=int(len(above)), names=list(above.sort_values("ps", ascending=False).province),
                           visit_share=float(100 * above.visits.sum() / V.sum()), receipt_share=float(100 * above.receipts.sum() / R.sum()))
    res["yield_extremes"] = dict(top=float(o.yield_thb[0]), bottom=float(o.yield_thb[76]), ratio=float(o.yield_thb[0] / o.yield_thb[76]),
                                 median_yield=float(df.yield_thb.median()), median_ps=float(df.ps.median()))

    # ---- Table 3: inequality with within-region bootstrap
    variables = [("Receipts (all)", "receipts"), ("Visits (all)", "visits"), ("Receipts, Thai visitors", "receipts_thai"),
                 ("Receipts, foreign visitors", "receipts_foreign"), ("Visits, Thai visitors", "visits_thai"), ("Visits, foreign visitors", "visits_foreign")]
    g = df.region.values
    idx_by_reg = [np.where(g == r)[0] for r in REGIONS]
    boot_idx = np.array([np.concatenate([rng.choice(ix, len(ix), replace=True) for ix in idx_by_reg]) for _ in range(a.boot)])
    g_b = np.concatenate([np.full(len(ix), k) for k, ix in enumerate(idx_by_reg)])
    rows, ineq = [], {}
    for lab, col in variables:
        x = df[col].values
        G, T, B = shape_stats(x, g)
        bs = np.array([shape_stats(x[bi], g_b) for bi in boot_idx])
        ci = np.percentile(bs, [2.5, 97.5], axis=0)
        ineq[col] = dict(gini=G, theil=T, between=B, top3=100 * topk(x, 3), top5=100 * topk(x, 5),
                         gini_ci=ci[:, 0].tolist(), theil_ci=ci[:, 1].tolist(), between_ci=ci[:, 2].tolist())
        rows.append([lab, f"{G:.3f}", f"[{ci[0, 0]:.3f}, {ci[1, 0]:.3f}]", f"{T:.3f}", f"[{ci[0, 1]:.3f}, {ci[1, 1]:.3f}]", f"{B:.1f}",
                     f"[{ci[0, 2]:.1f}, {ci[1, 2]:.1f}]", f"{100 * topk(x, 3):.1f}", f"{100 * topk(x, 5):.1f}"])
    res["inequality"] = ineq
    tex_table("table3", a.outdir, "lrcrcrcrr", r"Variable & Gini & Resampling interval & Theil $T$ & Resampling interval & Between (\%) & Resampling interval & Top 3 (\%) & Top 5 (\%)", rows)
    pc = df.receipts * 1e6 / df.population
    res["per_resident"] = dict(gini_equal=gini(pc.values), gini_popweighted=wgini(pc.values, df.population.values))

    # ---- decomposition (Eq. 8-9), Table 4
    def decomp(d):
        f = d.visits_foreign.sum() / d.visits.sum()
        yT = d.receipts_thai.sum() * 1e6 / d.visits_thai.sum()
        yF = d.receipts_foreign.sum() * 1e6 / d.visits_foreign.sum()
        y = d.receipts * 1e6 / d.visits
        fi = d.visits_foreign / d.visits
        yTi = d.receipts_thai * 1e6 / d.visits_thai
        yFi = d.receipts_foreign * 1e6 / d.visits_foreign
        comp = (fi - f) * (yF - yT)
        within = (1 - fi) * (yTi - yT) + fi * (yFi - yF)
        ybar_ = d.receipts.sum() * 1e6 / d.visits.sum()
        assert np.allclose(y - ybar_, comp + within, rtol=1e-9, atol=1e-6)               # identity of Eq. 8
        C1 = np.cov(y, comp, ddof=0)[0, 1] / np.var(y)
        C2 = np.cov(y, within, ddof=0)[0, 1] / np.var(y)
        return dict(comp_share=100 * C1, within_share=100 * C2, yT=yT, yF=yF, ybar=ybar_), y, comp, within, fi
    dec, y, comp, within, fi = decomp(df)
    dec_x, *_ = decomp(df[df.province != "Bangkok"])                              # benchmarks re-computed for the 76 provinces
    f77 = df.visits_foreign.sum() / df.visits.sum()                              # alternative: national (77-province) benchmarks kept fixed
    yT77, yF77 = dec["yT"], dec["yF"]
    d76 = df[df.province != "Bangkok"]
    y76 = d76.receipts * 1e6 / d76.visits
    comp76 = (d76.visits_foreign / d76.visits - f77) * (yF77 - yT77)
    share_fixed = 100 * np.cov(y76, comp76, ddof=0)[0, 1] / np.var(y76)
    res["decomposition"] = dict(**dec, comp_share_exBKK=dec_x["comp_share"], comp_share_exBKK_fixed_benchmarks=float(share_fixed), ratio_F_T=dec["yF"] / dec["yT"],
                                spearman_ps_foreign=float(stats.spearmanr(df.ps, df.f).statistic))
    sel = ["Phuket", "Krabi", "Bangkok", "Surat Thani", "Phangnga", "Chon Buri", "Chiang Mai", "Nong Bua Lam Phu"]
    rows = []
    for nme in sel:
        i = int(np.where(df.province == nme)[0][0])
        rows.append([nme, fmt(y[i]), f"{y[i] - dec['ybar']:,.0f}".replace("-", "$-$"), f"{comp[i]:,.0f}".replace("-", "$-$"),
                     f"{within[i]:,.0f}".replace("-", "$-$"), f"{100 * fi[i]:.1f}"])
    tex_table("table4", a.outdir, "lrrrrr", r"Province & Receipts/visit (THB) & Deviation from national (THB) & Composition effect (THB) & Within-type rate effect (THB) & Foreign visits (\%)", rows)

    # ---- spatial (Table 5, Figure 3)
    D = haversine(df.lat.values, df.lon.values)
    Dn = D.copy(); np.fill_diagonal(Dn, np.inf)
    res["spatial_setup"] = dict(min_band_no_islands_km=float(Dn.min(1).max()))
    band_km = float(Dn.min(1).max()) + 1e-9                    # smallest distance leaving no province without a neighbour
    res["spatial_setup"]["band_km"] = band_km
    specs = [("idw", "Inverse distance", None), ("knn", "5 nearest neighbours", None), ("band", "Distance band", band_km)]
    vars_ = [("ln yield", np.log(df.yield_thb.values)), ("Punch Score", df.ps.values), ("ln receipts", np.log(df.receipts.values)), ("Foreign visitor share", df.f.values)]
    rows, moran_res = [], []
    for vlab, x in vars_:
        for kind, wlab, band in specs:
            W = weights(D, kind, band)
            I, zsc, p = moran_perm(x, W, a.perm, rng)
            moran_res.append(dict(variable=vlab, weights=wlab, I=I, z=zsc, p=p))
            rows.append([vlab, wlab, f"{I:.3f}", f"{zsc:.2f}", "$<0.001$" if p < 0.001 else f"{p:.3f}"])
    res["moran"] = moran_res
    res["moran_null_expectation"] = -1.0 / 76
    tex_table("table5", a.outdir, "llrrr", r"Variable & Spatial weights & Moran's $I$ & $z$ & $p$", rows)
    W5 = weights(D, "knn")
    lny = np.log(df.yield_thb.values)
    quad, Ii, pl = local_moran(lny, W5, a.perm, rng)
    lab = {1: "High-High", 2: "Low-High", 3: "Low-Low", 4: "High-Low"}
    sig = pl < 0.05
    res["lisa"] = {lab[q]: list(df.province[(quad == q) & sig]) for q in (1, 2, 3, 4)}
    res["lisa"]["n_significant"] = int(sig.sum())
    res["lisa"]["n_not_significant"] = int((~sig).sum())
    freq = np.zeros(77)                                   # stability of the classification across independent permutation seeds
    n_rep = 20
    for r_ in range(n_rep):
        q_, _, p_ = local_moran(lny, W5, a.perm, np.random.default_rng(np.random.SeedSequence([a.seed, 777, r_])))
        freq += (p_ < 0.05)
    freq /= n_rep
    res["lisa"]["seed_frequency"] = {df.province[i]: float(freq[i]) for i in range(77) if sig[i] or freq[i] > 0}
    order = [i for q in (1, 2, 3, 4) for i in np.where((quad == q) & sig)[0]]
    near = [i for i in range(77) if (not sig[i]) and freq[i] >= 0.25]
    res["lisa"]["borderline"] = [df.province[i] for i in near]
    rows = [[df.province[i], df.region[i], lab[quad[i]], f"{Ii[i]:.2f}".replace("-", "$-$"), ("$<0.001$" if pl[i] < 0.001 else f"{pl[i]:.3f}"), f"{100 * freq[i]:.0f}"] for i in order + near]
    tex_table("table9", a.outdir, "llcrrr", r"Province & Region & Cluster type & Local $I_i$ & Pseudo-$p$ & Significant in (\% of 20 seeds)", rows)
    fig_sig, fig_quad = sig, quad

    # ---- Table 6: Spearman correlations with the Punch Score
    cvars = [("Foreign visitor share (\\%)", 100 * df.f), ("Average daily expenditure (THB/person/day)$^\\dagger$", df.expenditure_thb_day),
             ("Average length of stay (days)$^\\dagger$", df.los_days), ("Occupancy rate (\\%)", df.occupancy),
             ("Accommodation guests / visitors (\\%)", 100 * df.guests / df.visits), ("Rooms per 1,000 visits", 1000 * df.rooms / df.visits),
             ("GPP per capita (THB)", df.gpp_mthb * 1e6 / df.population), ("Visitors per resident", df.visits / df.population)]
    nb = (df.province != "Bangkok").values
    rows, corr = [], []
    pf = lambda p: "$<0.001$" if p < 0.001 else f"{p:.3f}"
    for lab_, x in cvars:
        r77 = stats.spearmanr(df.ps, x)
        r76 = stats.spearmanr(df.ps[nb], x[nb])
        corr.append(dict(variable=lab_, rho=r77.statistic, p=r77.pvalue, rho_exBKK=r76.statistic, p_exBKK=r76.pvalue))
        rows.append([lab_, f"{r77.statistic:.2f}", pf(r77.pvalue), f"{r76.statistic:.2f}", pf(r76.pvalue)])
    res["correlations"] = corr
    tex_table("table6", a.outdir, "lrrrr", r"Variable & $\rho$ (77) & $p$ & $\rho$ (excluding Bangkok) & $p$", rows)
    res["intensity"] = dict(national=100 * R.sum() / df.gpp_mthb.sum(),
                            national_published=100 * R.sum() / 18683885.0,
                            bangkok=float(100 * df.receipts[df.province == "Bangkok"].iloc[0] / df.gpp_mthb[df.province == "Bangkok"].iloc[0]),
                            phuket=float(100 * df.receipts[df.province == "Phuket"].iloc[0] / df.gpp_mthb[df.province == "Phuket"].iloc[0]),
                            krabi=float(100 * df.receipts[df.province == "Krabi"].iloc[0] / df.gpp_mthb[df.province == "Krabi"].iloc[0]))

    # ---- robustness excluding Bangkok
    d76 = df[nb].reset_index(drop=True)
    ps76 = punch_score(d76.receipts.values, d76.visits.values)
    g76 = d76.region.values
    gr, T_r, b_r = shape_stats(d76.receipts.values, g76)
    gv, T_v, b_v = shape_stats(d76.visits.values, g76)
    res["exBKK"] = dict(receipts=float(d76.receipts.sum()), visits=float(d76.visits.sum()), gini_r=gr, gini_v=gv, theil_r=T_r, theil_v=T_v,
                        between_r=b_r, between_v=b_v, n_above100=int((ps76 > 100).sum()), median_ps=float(np.median(ps76)))
    base = res["inequality"]
    res["baseline_row"] = dict(gini_r=base["receipts"]["gini"], gini_v=base["visits"]["gini"], theil_r=base["receipts"]["theil"], theil_v=base["visits"]["theil"],
                               between_r=base["receipts"]["between"], between_v=base["visits"]["between"], n_above100=res["above100"]["n"], median_ps=res["yield_extremes"]["median_ps"])
    b_, x_ = res["baseline_row"], res["exBKK"]
    rows = [["Total receipts (M THB)", fmt(R.sum(), 2), fmt(x_["receipts"], 2)],
            ["Total visits", fmt(V.sum()), fmt(x_["visits"])],
            ["Gini, receipts", f"{b_['gini_r']:.3f}", f"{x_['gini_r']:.3f}"],
            ["Gini, visits", f"{b_['gini_v']:.3f}", f"{x_['gini_v']:.3f}"],
            ["Theil $T$, receipts", f"{b_['theil_r']:.3f}", f"{x_['theil_r']:.3f}"],
            ["Theil $T$, visits", f"{b_['theil_v']:.3f}", f"{x_['theil_v']:.3f}"],
            ["Between-region share, receipts (\\%)", f"{b_['between_r']:.1f}", f"{x_['between_r']:.1f}"],
            ["Between-region share, visits (\\%)", f"{b_['between_v']:.1f}", f"{x_['between_v']:.1f}"],
            ["Provinces with Punch Score $>100$", b_["n_above100"], x_["n_above100"]],
            ["Median Punch Score", f"{b_['median_ps']:.1f}", f"{x_['median_ps']:.1f}"]]
    tex_table("table7", a.outdir, "lrr", r"Statistic & Baseline (77 provinces) & Excluding Bangkok (76 provinces)", rows)

    # ---- figures
    fig, ax = plt.subplots(figsize=(5.2, 5))
    for col, lab_, st in [("receipts", "Receipts", "-"), ("visits", "Recorded visits", "--")]:
        x = np.sort(df[col].values)
        L = np.concatenate([[0], np.cumsum(x) / x.sum()])
        ax.plot(np.linspace(0, 1, 78), L, st, label=f"{lab_} (Gini {gini(x):.3f})")
    ax.plot([0, 1], [0, 1], ":", color="grey", label="Line of equality")
    ax.set_xlabel("Cumulative share of provinces"); ax.set_ylabel("Cumulative share"); ax.legend(frameon=False); ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig/fig1_lorenz.pdf"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.2, 7))
    sc = ax.scatter(df.lon, df.lat, s=np.maximum(df.receipts / df.receipts.max() * 1400, 14), c=np.log10(df.ps / 100), cmap="RdBu_r", vmin=-0.8, vmax=0.8,
                    edgecolor="k", linewidth=0.3, alpha=0.9)
    for nme in ["Phuket", "Krabi", "Bangkok", "Surat Thani", "Chon Buri", "Chiang Mai", "Phangnga"]:
        r_ = df[df.province == nme].iloc[0]; ax.annotate(nme, (r_.lon, r_.lat), fontsize=7, xytext=(4, 4), textcoords="offset points")
    plt.colorbar(sc, label=r"$\log_{10}(PS/100)$", shrink=0.6); ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig/fig2_yield_map.pdf"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.2, 7))
    colors = {1: "#d7191c", 2: "#fdae61", 3: "#2c7bb6", 4: "#abd9e9"}
    ax.scatter(df.lon[~sig], df.lat[~sig], s=22, c="lightgrey", edgecolor="grey", linewidth=0.3, label="Not significant")
    for q in (1, 2, 3, 4):
        m_ = sig & (quad == q)
        if m_.any():
            ax.scatter(df.lon[m_], df.lat[m_], s=60, c=colors[q], edgecolor="k", linewidth=0.4, label=lab[q])
            for _, r_ in df[m_].iterrows():
                ax.annotate(r_.province, (r_.lon, r_.lat), fontsize=6.5, xytext=(4, 3), textcoords="offset points")
    ax.legend(frameon=False, fontsize=8); ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig/fig3_lisa.pdf"); plt.close(fig)

    json.dump(res, open(f"{a.outdir}/results.json", "w"), indent=1, default=float)
    man = dict(seed=a.seed, boot=a.boot, perm=a.perm, data_sha256=hashlib.sha256(open(a.data, "rb").read()).hexdigest(),
               code_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest(), python=sys.version.split()[0], numpy=np.__version__,
               scipy=scipy.__version__, pandas=pd.__version__, platform=platform.platform(), runtime_s=round(time.time() - t0, 1))
    man["run_id"] = hashlib.sha256((man["data_sha256"] + man["code_sha256"] + f"{a.seed}{a.boot}{a.perm}").encode()).hexdigest()[:12]
    json.dump(man, open(f"{a.outdir}/manifest.json", "w"), indent=1)
    print("done", man["run_id"], man["runtime_s"], "s")


if __name__ == "__main__":
    main()
