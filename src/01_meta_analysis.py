#!/usr/bin/env python3
"""Primary meta-analysis: per-study estimates, pooling, power, FPRP, FST.

Produces Supplementary Tables S2 (derived counts), S3 (estimates),
S4 (power and FPRP). Run from the repository root:  python src/01_meta_analysis.py
"""
import sys
from math import exp
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eaat.data import load_genotypes
from eaat.stats import (alleles, fmt_or, fprp, fst_two_pop, hwe, meta,
                        n_required, or_2x2, power_allelic)

G = load_genotypes()
MDD = ["rodek22", "thaweethee19"]          # eligible for pooling
RISK = "C"                                  # low-promoter-activity allele


def models(cs, ct):
    """Allelic, recessive and dominant log-OR + SE."""
    ac, at = alleles(cs), alleles(ct)
    return {
        "allelic":   or_2x2(ac["C"], ac["A"], at["C"], at["A"]),
        "recessive": or_2x2(cs["CC"], cs["AA"] + cs["AC"],
                            ct["CC"], ct["AA"] + ct["AC"]),
        "dominant":  or_2x2(cs["CC"] + cs["AC"], cs["AA"],
                            ct["CC"] + ct["AC"], ct["AA"]),
    }


print("=" * 78)
print("PER-STUDY ESTIMATES   (risk allele = C)")
print("=" * 78)
for sid, rec in G.items():
    cs, ct = rec["case"], rec["control"]
    ac, at = alleles(cs), alleles(ct)
    print(f"\n{rec['meta']['study']} {rec['meta']['year']} "
          f"({rec['meta']['phenotype']}, {rec['meta']['ancestry']})"
          f"{'   [shared controls]' if rec['meta']['shared_controls'] else ''}")
    print(f"  n cases={sum(cs.values())}  n controls={sum(ct.values())}")
    print(f"  C freq  cases={ac['C'] / (2 * sum(cs.values())):.3f}  "
          f"controls={at['C'] / (2 * sum(ct.values())):.3f}")
    for grp, g in (("cases", cs), ("controls", ct)):
        chi, p = hwe(g)
        print(f"  HWE {grp:<9} chi2={chi:.2f}  p={p:.3f}")
    for name, (lo, se) in models(cs, ct).items():
        print(f"  {name:<10} OR = {fmt_or(lo, se)}")

print("\n" + "=" * 78)
print("POOLED ESTIMATES  (MDD only; the SZ series shares its controls)")
print("=" * 78)
for name in ("allelic", "recessive", "dominant"):
    ys, ses = [], []
    for sid in MDD:
        lo, se = models(G[sid]["case"], G[sid]["control"])[name]
        ys.append(lo); ses.append(se)
    r = meta(ys, ses)
    from scipy import stats as st
    pf = 2 * (1 - st.norm.cdf(abs(r["mu_f"] / r["se_f"])))
    lo_hk = r["mu_r"] - r["t"] * r["se_hk"]
    hi_hk = r["mu_r"] + r["t"] * r["se_hk"]
    print(f"\n{name.upper()}  (k = {r['k']})")
    print(f"  Fixed  (IV)    OR = {fmt_or(r['mu_f'], r['se_f'], 3)}  p = {pf:.4f}")
    print(f"  Random (DL)    OR = {fmt_or(r['mu_r'], r['se_r'], 3)}")
    print(f"  Random (HK t)  OR = {exp(r['mu_r']):.3f} "
          f"({exp(lo_hk):.3f}-{exp(hi_hk):.3f})   <-- report this")
    print(f"  Q = {r['Q']:.2f}  p = {r['pQ']:.4f}  "
          f"I2 = {r['I2']:.1f}%  tau2 = {r['tau2']:.4f}")

print("\n" + "=" * 78)
print("POWER")
print("=" * 78)
print(f"{'true OR':>9} | " + " | ".join(
    f"{G[s]['meta']['study'].split()[0]:>12}" for s in G))
for OR in (1.02, 1.05, 1.10, 1.20, 1.57):
    cells = []
    for sid, rec in G.items():
        at = alleles(rec["control"])
        maf = at["C"] / (2 * sum(rec["control"].values()))
        cells.append(f"{power_allelic(sum(rec['case'].values()),
                     sum(rec['control'].values()), maf, OR) * 100:>11.1f}%")
    print(f"{OR:>9.2f} | " + " | ".join(cells))

print("\nCases needed for 80% power at MAF 0.41:")
for OR in (1.05, 1.10, 1.20):
    print(f"  OR {OR:.2f}  alpha=0.05 : {n_required(0.41, OR, 0.05):>10,}"
          f"   alpha=5e-8 : {n_required(0.41, OR, 5e-8):>10,}")

print("\n" + "=" * 78)
print("FALSE-POSITIVE REPORT PROBABILITY  (recessive finding, p ~ 0.001)")
print("power evaluated at the PLAUSIBLE true OR, not the observed one")
print("=" * 78)
rodek = G["rodek22"]
maf = alleles(rodek["control"])["C"] / (2 * sum(rodek["control"].values()))
print(f"{'prior':>8} | " + " | ".join(f"{f'OR {o}':>9}" for o in (1.05, 1.20, 2.38)))
for prior in (0.25, 0.05, 0.01, 0.001):
    cells = [f"{fprp(0.001, power_allelic(161, 462, maf, o, alpha=0.001), prior):>9.3f}"
             for o in (1.05, 1.20, 2.38)]
    print(f"{prior:>8.3f} | " + " | ".join(cells))

print("\n" + "=" * 78)
print("ANCESTRY DIVERGENCE")
print("=" * 78)
f_pl = alleles(G["rodek22"]["control"]); f_th = alleles(G["thaweethee19"]["control"])
p1 = f_pl["C"] / (f_pl["C"] + f_pl["A"]); p2 = f_th["C"] / (f_th["C"] + f_th["A"])
print(f"  C freq Poland = {p1:.3f}   C freq Thailand = {p2:.3f}")
print(f"  two-population Wright's FST = {fst_two_pop(p1, p2):.4f}")
