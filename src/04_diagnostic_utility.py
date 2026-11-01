#!/usr/bin/env python3
"""Ceiling on diagnostic discrimination and predictive value.

Produces Supplementary Table S6.  Run:  python src/04_diagnostic_utility.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eaat.stats import auc_additive

RISK_FREQ, PREVALENCE = 0.41, 0.20

print("=" * 72)
print("DISCRIMINATION: AUC of rs4354668 alone (log-additive, freq 0.41)")
print("=" * 72)
print(f"{'per-allele OR':>14} | {'AUC':>7} | interpretation")
for OR in (1.02, 1.05, 1.10, 1.20, 1.57, 2.00, 3.00):
    a = auc_additive(RISK_FREQ, OR)
    tag = ("indistinguishable from chance" if a < 0.55 else
           "no clinical value" if a < 0.60 else
           "poor" if a < 0.70 else "moderate")
    mark = "  <-- published estimate" if OR == 1.57 else ""
    print(f"{OR:>14.2f} | {a:>7.3f} | {tag}{mark}")
print("\n  Chance = 0.500. Screening instruments are unusable below ~0.75.")
print("  The entire depression polygenic score reaches only ~0.57-0.60.")

print("\n" + "=" * 72)
print("PREDICTIVE VALUE of the CC genotype, taking the OR at face value")
print(f"Lifetime prevalence {PREVALENCE}; CC frequency {RISK_FREQ**2:.2f}")
print("=" * 72)
print(f"{'recessive OR':>13} | {'risk if CC':>11} | {'if not CC':>10} | "
      f"{'sens':>6} | {'spec':>6}")
f_hom = RISK_FREQ ** 2
for OR in (2.38, 1.50, 1.20, 1.05):
    r0 = PREVALENCE / (f_hom * OR + (1 - f_hom))
    r_pos, r_neg = r0 * OR, r0
    sens = f_hom * r_pos / PREVALENCE
    spec = 1 - (f_hom * (1 - r_pos)) / (1 - PREVALENCE)
    print(f"{OR:>13.2f} | {r_pos*100:>10.1f}% | {r_neg*100:>9.1f}% | "
          f"{sens*100:>5.1f}% | {spec*100:>5.1f}%")
print("\n  Even at the published OR, genotyping shifts lifetime risk from ~16%")
print("  to ~39% while missing ~68% of future cases. A single question about")
print("  family history outperforms it.")
