#!/usr/bin/env python3
"""Replication ledger, sign test, excess significance, fail-safe N.

Produces Supplementary Table S5.  Run:  python src/03_ledger_bias.py
"""
import sys
from math import log, sqrt
from pathlib import Path

from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eaat.data import load_ledger
from eaat.stats import power_allelic

L = load_ledger()
PGC_CASES = 688_808

print("=" * 80)
print("REPLICATION LEDGER -- every located case-control test of rs4354668")
print("=" * 80)
print(f"{'study':<20}|{'dx':<5}|{'ancestry':<9}|{'cases':>6}|{'ctrls':>6}|"
      f"{'sig':>4}|{'counts':>7}")
for r in L:
    approx = "~" if r["n_approximate"] else " "
    print(f"{r['study'] + ' ' + str(r['year']):<20}|{r['phenotype']:<5}|"
          f"{r['ancestry']:<9}|{approx}{r['n_cases']:>5}|{approx}{r['n_controls']:>5}|"
          f"{'YES' if r['significant'] else 'no':>4}|"
          f"{'yes' if r['counts_published'] else 'NO':>7}")
tot_ca = sum(r["n_cases"] for r in L)
mdd_ca = sum(r["n_cases"] for r in L if r["phenotype"] == "MDD")
print(f"\n  All psychiatric phenotypes: {tot_ca:,} cases / "
      f"{sum(r['n_controls'] for r in L):,} controls")
print(f"  MDD only:                   {mdd_ca:,} cases / "
      f"{sum(r['n_controls'] for r in L if r['phenotype']=='MDD'):,} controls")
print(f"  The entire rs4354668 MDD literature is 1/{PGC_CASES // mdd_ca:,}th "
      f"of one PGC GWAS.")

print("\n" + "=" * 80)
print("TEST 1 -- SIGN TEST ON DIRECTION")
print("=" * 80)
dirs = [r["direction"] for r in L if r["direction"] != 0]
k_pos = sum(1 for d in dirs if d > 0)
p_sign = stats.binomtest(k_pos, len(dirs), 0.5, alternative="two-sided").pvalue
print(f"  {k_pos}/{len(dirs)} same direction; exact binomial p = {p_sign:.3f}")
print("  At k=4 this cannot reach significance even at 4/4, and a literature")
print("  reporting only confirmatory analyses produces this by construction.")

print("\n" + "=" * 80)
print("TEST 2 -- EXCESS SIGNIFICANCE (Ioannidis & Trikalinos 2007)")
print("=" * 80)
O = sum(1 for r in L if r["significant"])
print(f"  Observed significant findings O = {O} of {len(L)}\n")
print(f"{'true OR':>9} | {'expected E':>11} | {'O-E':>6} | {'binom p':>9} | verdict")
for OR in (1.02, 1.05, 1.10, 1.20, 1.45):
    E = sum(power_allelic(r["n_cases"], r["n_controls"],
                          r["control_freq_C"], OR) for r in L)
    p = stats.binomtest(O, len(L), E / len(L), alternative="greater").pvalue
    v = "EXCESS" if p < 0.10 else "consistent" if p < 0.5 else "under-reported"
    print(f"{OR:>9.2f} | {E:>11.2f} | {O-E:>+6.2f} | {p:>9.3f} | {v}")
print("\n  At the effect sizes modern psychiatric genetics observes (OR 1.02-1.10)")
print("  these six studies had power to expect well under one significant")
print("  result. Two were published. The excess resolves only by assuming the")
print("  conclusion under test.")

print("\n" + "=" * 80)
print("TEST 3 -- FAIL-SAFE N against an equivalence bound of OR 1.05")
print("=" * 80)
mu, se, target = log(1.449), 0.1116, log(1.05)
for n in (200, 500, 1000):
    se_k = sqrt(4 / (2 * n * 0.41 * 0.59))
    for kk in range(1, 400):
        if ((1/se**2) * mu) / ((1/se**2) + kk / se_k**2) < target:
            print(f"  Null studies of {n:>4} cases each needed: {kk}")
            break
print("  The pooled estimate rests on 261 cases. A handful of null studies")
print("  erases it; the PGC data hold 2,600x that, unanalysed for this SNP.")
