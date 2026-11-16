#!/usr/bin/env python3
"""
build_dataset.py -- extracts the 77-province 2024 analysis dataset from the compiled data tables (Tables 1-3 of the data
PDF, which reproduce the Ministry of Tourism and Sports revised workbook of 28 August 2025, the Domestic Tourism
Statistics workbooks, NESDC GPP, NSO population/area, and OCHA-derived centroids) and validates it against the
published regional and national totals.

    python build_dataset.py data_tourism_2024.pdf provinces_2024.csv        (needs pdftotext)
"""
import csv
import re
import subprocess
import sys

pdf, out = sys.argv[1], sys.argv[2]
txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True, check=True).stdout
T2, T3 = "\u0e15\u0e32\u0e23\u0e32\u0e07\u0e17\u0e35\u0e48 2", "\u0e15\u0e32\u0e23\u0e32\u0e07\u0e17\u0e35\u0e48 3"   # Thai headings "Table 2", "Table 3"
t1, rest = txt.split(T2, 1)
t2, t3 = rest.split(T3, 1)
num = r"([\d,]+(?:\.\d+)?)"
name = r"^\s*([A-Za-z][A-Za-z .()\-]*?(?:\(\d+\))?)\s+"
f = lambda s: float(s.replace(",", ""))
R1 = re.compile(name + r"\s+".join([num] * 8) + r"\s*$")
R2 = re.compile(name + r"\s+".join([num] * 6) + r"\s*$")
R3 = re.compile(name + r"\s+".join([r"([\d,]+\.\d+)"] * 4 + [r"(?:(\d+\.\d+)|\u2013)"] * 2) + r"\s*$")


def rows(text, rx):
    out = []
    for line in text.splitlines():
        m = rx.match(line)
        if m:
            out.append((m.group(1).strip(), [g for g in m.groups()[1:]]))
    return out


a, b, c = rows(t1, R1), rows(t2, R2), rows(t3, R3)
REGIONS = {"Central (incl. West) (17)": "Central", "East (8)": "East", "North (17)": "North", "Northeast (20)": "Northeast", "South (14)": "South"}
prov, region, cur = [], {}, "Bangkok"
sub = {}
for (n1, v1), (n2, v2), (n3, v3) in zip(a, b, c):
    assert n1 == n2 == n3, (n1, n2, n3)
    if n1 in REGIONS or n1 == "Thailand (77)":
        sub[n1] = (v1, v2, v3)
        cur = REGIONS.get(n1, cur)
        continue
    region[n1] = cur
    prov.append((n1, cur, [f(x) for x in v1], [f(x) for x in v2], [None if x is None else f(x) for x in v3]))
assert len(prov) == 77 and len({p[0] for p in prov}) == 77, len(prov)
counts = {r: sum(1 for p in prov if p[1] == r) for r in ["Bangkok", "Central", "East", "North", "Northeast", "South"]}
assert counts == {"Bangkok": 1, "Central": 17, "East": 8, "North": 17, "Northeast": 20, "South": 14}, counts
# ---- validation against the published regional / national rows
cols1 = ["visits", "visits_thai", "visits_foreign", "receipts", "receipts_thai", "receipts_foreign", "guests", "occupancy"]
ok = True
for reg, label in list(REGIONS.items()) + [("Thailand (77)", None)]:
    members = [p for p in prov if (label is None or p[1] == label)]
    for k, nme in enumerate(cols1[:7]):
        s = sum(p[2][k] for p in members)
        pub = f(sub[reg][0][k])
        if abs(s - pub) > (0.02 if "receipts" in nme else 0.5):
            ok = False
            print("MISMATCH", reg, nme, s, pub)
    rooms = sum(p[3][0] for p in members); pub_r = f(sub[reg][1][0])
    if abs(rooms - pub_r) > 0.5:
        ok = False; print("MISMATCH rooms", reg, rooms, pub_r)
print("77 provinces; regional and national totals reconcile:", ok)
bad = [p[0] for p in prov if abs(p[2][1] + p[2][2] - p[2][0]) > 0.5 or abs(p[2][4] + p[2][5] - p[2][3]) > 0.02]
print("Thai+foreign != total in:", bad)
gpp_sum = sum(p[3][3] for p in prov)
print("GPP: provincial sum %.0f vs published national %.0f" % (gpp_sum, f(sub["Thailand (77)"][1][3])))
hdr = ["province", "region"] + cols1 + ["rooms", "los_days", "expenditure_thb_day", "gpp_mthb", "population", "area_km2",
                                         "los_thai", "los_foreign", "exp_thai", "exp_foreign", "lat", "lon"]
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(hdr)
    for n, r, v1, v2, v3 in prov:
        w.writerow([n, r] + v1 + v2 + v3)
print("wrote", out)
