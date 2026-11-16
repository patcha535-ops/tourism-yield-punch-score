#!/usr/bin/env python3
"""
verify_claims.py -- compares every number quoted in the submitted manuscript (hard-coded below, as printed) with the values
recomputed by te_analysis.py.  A number "matches" if it agrees after rounding to the precision printed in the paper.

    python verify_claims.py results_te/results.json
"""
import json
import sys

import pandas as pd

r = json.load(open(sys.argv[1]))
df = pd.read_csv("provinces_2024.csv")
rows = []


def chk(label, paper, mine, dec=None, tol=None):
    ok = abs(paper - mine) <= (tol if tol is not None else 0.5 * 10 ** (-dec) + 1e-9)
    rows.append((ok, label, paper, mine))


reg = r["regions"]
chk("Bangkok+South visit share (%)", 33.2, reg["Bangkok"]["vshare"] + reg["South"]["vshare"], tol=0.051)
chk("Bangkok+South receipts share (%)", 65.9, reg["Bangkok"]["rshare"] + reg["South"]["rshare"], tol=0.051)
chk("Central+Northeast visit share (%)", 40.2, reg["Central"]["vshare"] + reg["Northeast"]["vshare"], tol=0.051)
chk("Central+Northeast receipts share (%)", 11.6, reg["Central"]["rshare"] + reg["Northeast"]["rshare"], tol=0.051)
for k, v in dict(South=203.1, Bangkok=193.7, East=100.0, North=65.1, Central=29.5, Northeast=27.9).items():
    chk(f"Regional Punch Score {k}", v, reg[k]["ps"], 1)
chk("National receipts/visit (THB)", 8272, r["national"]["ybar"], 0)
chk("Foreign share of visits (%)", 23.5, 100 * r["national"]["foreign_visit_share"], 1)
chk("Foreign share of receipts (%)", 61.5, 100 * r["national"]["foreign_receipt_share"], 1)
chk("Receipts per Thai visit", 4162, r["national"]["ybar_thai"], 0)
chk("Receipts per foreign visit", 21632, r["national"]["ybar_foreign"], 0)
chk("Foreign/Thai yield ratio", 5.2, r["decomposition"]["ratio_F_T"], 1)
ab = r["above100"]
chk("Provinces above benchmark", 7, ab["n"], 0)
chk("Their share of visits (%)", 37.5, ab["visit_share"], 1)
chk("Their share of receipts (%)", 77.2, ab["receipt_share"], 1)
ex = r["yield_extremes"]
chk("Highest receipts/visit (Phuket)", 39387, ex["top"], 0); chk("Lowest receipts/visit", 1172, ex["bottom"], 0)
chk("Highest/lowest ratio", 34, ex["ratio"], 0); chk("Median yield (THB)", 2416, ex["median_yield"], 0); chk("Median Punch Score", 29.2, ex["median_ps"], 1)
for nme, ps in dict(Phuket=476.2, Krabi=197.3, Bangkok=193.7, **{"Nong Bua Lam Phu": 14.2}).items():
    i = r["ranking"]["province"].index(nme); chk(f"Punch Score {nme}", ps, r["ranking"]["ps"][i], 1)
I = r["inequality"]
paper = {"receipts": (0.827, 1.765, 60.0, 63.5, 71.7), "visits": (0.558, 0.639, 54.3, 28.8, 36.1), "receipts_thai": (0.668, 0.992, 58.5, 39.0, 46.6),
         "receipts_foreign": (0.938, 2.562, 60.8, 80.2, 89.6), "visits_thai": (0.502, 0.474, 49.6, 23.1, 30.9), "visits_foreign": (0.865, 1.839, 59.7, 59.9, 70.1)}
for k, (g, t, b, t3, t5) in paper.items():
    chk(f"Gini {k}", g, I[k]["gini"], 3); chk(f"Theil {k}", t, I[k]["theil"], 3); chk(f"Between % {k}", b, I[k]["between"], 1)
    chk(f"Top-3 % {k}", t3, I[k]["top3"], 1); chk(f"Top-5 % {k}", t5, I[k]["top5"], 1)
paper_ci = {"receipts": ((0.770, 0.851), (1.506, 2.008), (54.7, 85.1)), "visits": ((0.501, 0.601), (0.536, 0.745), (48.6, 71.1))}
for k, (gc, tc, bc) in paper_ci.items():                       # bootstrap intervals: Monte Carlo agreement within 0.01 / 0.05 / 1.5
    chk(f"Gini CI low {k}", gc[0], I[k]["gini_ci"][0], tol=0.01); chk(f"Gini CI high {k}", gc[1], I[k]["gini_ci"][1], tol=0.01)
    chk(f"Theil CI low {k}", tc[0], I[k]["theil_ci"][0], tol=0.05); chk(f"Theil CI high {k}", tc[1], I[k]["theil_ci"][1], tol=0.05)
    chk(f"Between CI low {k}", bc[0], I[k]["between_ci"][0], tol=1.5); chk(f"Between CI high {k}", bc[1], I[k]["between_ci"][1], tol=1.5)
chk("Gini receipts per resident, equal weights", 0.787, r["per_resident"]["gini_equal"], 3)
chk("Gini receipts per resident, population-weighted", 0.745, r["per_resident"]["gini_popweighted"], 3)
D = r["decomposition"]
chk("Composition share (%)", 44.6, D["comp_share"], 1); chk("Within-type share (%)", 55.4, D["within_share"], 1)
chk("Composition share excl. Bangkok, national benchmarks fixed (%)", 44.2, D["comp_share_exBKK_fixed_benchmarks"], 1)
chk("Spearman PS vs foreign share", 0.62, D["spearman_ps_foreign"], 2)
pm = {("ln yield", "Inverse distance"): (0.168, 10.56), ("ln yield", "5 nearest neighbours"): (0.434, 7.09), ("ln yield", "Distance band"): (0.494, 5.74),
      ("Punch Score", "Inverse distance"): (0.112, 8.84), ("Punch Score", "5 nearest neighbours"): (0.303, 6.10), ("Punch Score", "Distance band"): (0.415, 5.99),
      ("ln receipts", "Inverse distance"): (0.095, 6.33), ("ln receipts", "5 nearest neighbours"): (0.256, 4.15), ("ln receipts", "Distance band"): (0.273, 3.24),
      ("Foreign visitor share", "Inverse distance"): (0.175, 10.98), ("Foreign visitor share", "5 nearest neighbours"): (0.373, 6.12), ("Foreign visitor share", "Distance band"): (0.537, 6.32)}
for m in r["moran"]:
    pI, pz = pm[(m["variable"], m["weights"])]
    chk(f"Moran I {m['variable']} / {m['weights']}", pI, m["I"], 3)
    chk(f"Moran z {m['variable']} / {m['weights']} (permutation noise, tol 0.5)", pz, m["z"], tol=0.5)
chk("Distance band (km), reported 112", 112, r["spatial_setup"]["band_km"], tol=0.5)
L = r["lisa"]
chk("Number of High-High clusters", 8, len(L["High-High"]), 0); chk("Number of Low-Low clusters", 2, len(L["Low-Low"]), 0)
chk("Low-High provinces (Phatthalung)", 1, len(L["Low-High"]), 0)
corr = {"Foreign visitor share": (0.62, 0.61), "Average daily expenditure": (0.97, 0.97), "Average length of stay": (0.72, 0.71), "Occupancy rate": (0.60, 0.59),
        "Accommodation guests": (0.77, 0.77), "Rooms per 1,000": (0.73, 0.73), "GPP per capita": (0.27, 0.24), "Visitors per resident": (0.41, 0.39)}
for c in r["correlations"]:
    key = next(k for k in corr if c["variable"].startswith(k)); chk(f"Spearman {key}", corr[key][0], c["rho"], 2); chk(f"Spearman {key} excl. BKK", corr[key][1], c["rho_exBKK"], 2)
chk("Receipts / GPP national (%)", 15.7, r["intensity"]["national_published"], 1); chk("Receipts / GPP Bangkok (%)", 14.8, r["intensity"]["bangkok"], 1)
chk("Receipts / GPP Phuket (%)", 210, r["intensity"]["phuket"], 0); chk("Receipts / GPP Krabi (%)", 111, r["intensity"]["krabi"], 0)
x = r["exBKK"]
for lab, pv, mv, dd in [("Gini receipts excl. BKK", 0.773, x["gini_r"], 3), ("Gini visits excl. BKK", 0.492, x["gini_v"], 3), ("Theil receipts excl. BKK", 1.454, x["theil_r"], 3),
                        ("Theil visits excl. BKK", 0.427, x["theil_v"], 3), ("Between receipts excl. BKK", 28.3, x["between_r"], 1), ("Between visits excl. BKK", 17.9, x["between_v"], 1),
                        ("Provinces >100 excl. BKK", 9, x["n_above100"], 0), ("Median PS excl. BKK", 35.7, x["median_ps"], 1), ("Receipts excl. BKK", 1982313.45, x["receipts"], 2), ("Visits excl. BKK", 294752564, x["visits"], 0)]:
    chk(lab, pv, mv, dd)
tab4 = {"Phuket": (39387, 31115, 9212, 21903), "Krabi": (16316, 8045, 5626, 2419), "Bangkok": (16023, 7751, 3939, 3811), "Surat Thani": (14543, 6272, 4189, 2083),
        "Phangnga": (13033, 4761, 6618, -1857), "Chon Buri": (12306, 4034, 3210, 824), "Chiang Mai": (9312, 1041, 1287, -246), "Nong Bua Lam Phu": (1172, -7100, -3951, -3148)}
ybar = r["decomposition"]["ybar"]; fT, yT, yF = r["national"]["foreign_visit_share"], r["decomposition"]["yT"], r["decomposition"]["yF"]
for nme, (y_, dv, cp, wt) in tab4.items():
    d = df[df.province == nme].iloc[0]
    y = d.receipts * 1e6 / d.visits; fi = d.visits_foreign / d.visits
    comp = (fi - fT) * (yF - yT); within = (1 - fi) * (d.receipts_thai * 1e6 / d.visits_thai - yT) + fi * (d.receipts_foreign * 1e6 / d.visits_foreign - yF)
    chk(f"Table 4 {nme} yield", y_, y, 0); chk(f"Table 4 {nme} deviation", dv, y - ybar, 0); chk(f"Table 4 {nme} composition", cp, comp, 0); chk(f"Table 4 {nme} within", wt, within, 0)
# new statements introduced in the revision (all taken from the recomputation)
chk("Distance band threshold (km), revised text 111.6", 111.6, r["spatial_setup"]["band_km"], 1)
chk("Composition share excl. Bangkok, benchmarks recomputed (%), revised text 41.9", 41.9, D["comp_share_exBKK"], 1)
chk("Local clusters significant (revised text: 11)", 11, L["n_significant"], 0)
chk("Provinces without significant local pattern (revised text: 66)", 66, L["n_not_significant"], 0)
chk("Khon Kaen significant in 65% of permutation runs (revised text)", 0.65, L["seed_frequency"]["Khon Kaen"], 2)
bad = [x_ for x_ in rows if not x_[0]]
print(f"{len(rows) - len(bad)} of {len(rows)} quoted numbers reproduce")
for ok, lab, p, m in bad:
    print(f"  DIFFERS  {lab}: paper {p}  recomputed {m:.4f}")
