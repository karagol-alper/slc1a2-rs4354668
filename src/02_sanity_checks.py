#!/usr/bin/env python3
"""Six forensic checks on the positive finding.

Produces Supplementary Table S3 (quality control, case-case, sex interaction)
and the multiplicity accounting. Run:  python src/02_sanity_checks.py
"""
import sys
from math import sqrt
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eaat.data import SHARED_CONTROLS, load_genotypes
from eaat.stats import alleles, fmt_or, hwe, or_2x2, wright_f

G = load_genotypes()
mdd, sz = G["rodek22"]["case"], G["merk19"]["case"]
ctl = G["rodek22"]["control"]

print("=" * 78)
print("CHECK 1 -- ARE THE TWO EUROPEAN STUDIES INDEPENDENT?")
print("=" * 78)
same = G["rodek22"]["control"] == G["merk19"]["control"]
print(f"  Rodek 2022 controls : {G['rodek22']['control']}")
print(f"  Merk  2019 controls : {G['merk19']['control']}")
print(f"  Identical: {same}   shared-control set: {sorted(SHARED_CONTROLS)}")
print("  Same centre, same bioethics approval, same primers, same BcnI digest.")
print("  => ONE control series. Pooling them double-counts 462 people.")

print("\n" + "=" * 78)
print("CHECK 2 -- GENOTYPING QC: HETEROZYGOTE DEFICIT")
print("=" * 78)
print(f"{'series':<30} | {'n':>4} | {'HWE p':>6} | {'F':>7} | {'95% CI':>18}")
series = [("Rodek MDD cases (PL)", mdd), ("Merk SZ cases (PL)", sz),
          ("Rodek/Merk controls (PL)", ctl),
          ("Thaweethee MDD cases (TH)", G["thaweethee19"]["case"]),
          ("Thaweethee controls (TH)", G["thaweethee19"]["control"])]
for name, g in series:
    _, p = hwe(g)
    f, (lo, hi) = wright_f(g)
    flag = "  <-- largest" if f > 0.10 else ""
    print(f"{name:<30} | {sum(g.values()):>4} | {p:>6.3f} | {f:>+7.3f} | "
          f"[{lo:+.3f},{hi:+.3f}]{flag}")
print("\n  All pass the conventional HWE screen. The one series producing a")
print("  positive result carries the largest heterozygote deficit. The CI")
print("  includes zero: hypothesis-generating, not probative.")

print("\n" + "=" * 78)
print("CHECK 3 -- SAME LAB, SAME CONTROLS, 2.4x MORE PATIENTS")
print("=" * 78)
for lab, cases in (("MDD (n=161)", mdd), ("SZ  (n=393)", sz)):
    ac, at = alleles(cases), alleles(ctl)
    la, sa = or_2x2(ac["C"], ac["A"], at["C"], at["A"])
    lr, sr = or_2x2(cases["CC"], cases["AA"] + cases["AC"],
                    ctl["CC"], ctl["AA"] + ctl["AC"])
    p = 2 * (1 - stats.norm.cdf(abs(lr / sr)))
    print(f"  {lab}: CC={cases['CC']/sum(cases.values())*100:4.1f}%  "
          f"allelic OR={fmt_or(la,sa)}  recessive OR={fmt_or(lr,sr)}  p={p:.4f}")
print(f"  Controls:    CC={ctl['CC']/sum(ctl.values())*100:4.1f}%")
print("  The better-powered SZ series sits on the controls. The controls are")
print("  not the anomaly; the 161-patient MDD series is.")

print("\n" + "=" * 78)
print("CHECK 4 -- CASE-CASE TEST (removes the controls entirely)")
print("=" * 78)
lo, se = or_2x2(mdd["CC"], mdd["AA"] + mdd["AC"], sz["CC"], sz["AA"] + sz["AC"])
print(f"  CC genotype, MDD vs SZ: OR = {fmt_or(lo,se)}  "
      f"p = {2*(1-stats.norm.cdf(abs(lo/se))):.4f}")

print("\n" + "=" * 78)
print("CHECK 5 -- IS THE REPORTED FEMALE-SPECIFIC EFFECT REAL?")
print("Rodek compared p-values across strata. The test compares the ORs.")
print("=" * 78)
lf, sf = or_2x2(35, 81, 34, 204)      # women,  recessive, Rodek Table 5
lm, sm = or_2x2(14, 31, 40, 184)      # men
d, sd = lf - lm, sqrt(sf**2 + sm**2)
print(f"  Women: OR = {fmt_or(lf,sf)}")
print(f"  Men:   OR = {fmt_or(lm,sm)}")
print(f"  Ratio of ORs = {fmt_or(d,sd)}")
print(f"  Interaction: z = {d/sd:.2f}, p = {2*(1-stats.norm.cdf(abs(d/sd))):.2f}")
print("  No evidence of sex modification (Gelman & Stern 2006).")

print("\n" + "=" * 78)
print("CHECK 6 -- MULTIPLICITY")
print("=" * 78)
tbl = np.array([[mdd["CC"], mdd["AA"] + mdd["AC"]],
                [ctl["CC"], ctl["AA"] + ctl["AC"]]])
_, p_rec, _, _ = stats.chi2_contingency(tbl, correction=False)
print(f"  Unadjusted recessive p = {p_rec:.2e}")
for k, lab in ((45, "3 SNPs x 5 models x 3 strata"),
               (69, "+ 8 haplotypes x 3 strata"),
               (105, "+ 4 clinical vars x 3 effects x 3 SNPs")):
    print(f"  Bonferroni over {k:>3} [{lab}]: p = {min(1.0, p_rec*k):.3f}")
print("  Survives the grid. Does NOT survive the selection step: the model")
print("  was chosen by BIC, then inference reported from the chosen model.")
