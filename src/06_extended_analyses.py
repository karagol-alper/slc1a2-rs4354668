#!/usr/bin/env python3
"""Extended analyses: exact tests, Bayesian reanalysis, quantitative bias
analysis, vibration of effects, and a simulation of the literature under a
true null with selective reporting.

Produces Supplementary Tables S7-S11.  Run:  python src/06_extended_analyses.py
"""
import sys
from itertools import product
from math import exp, log, sqrt
from pathlib import Path

import numpy as np
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eaat.data import load_genotypes, load_ledger
from eaat.stats import alleles, meta, or_2x2, power_allelic

G = load_genotypes()
L = load_ledger()
RODEK_C, RODEK_T = G["rodek22"]["case"], G["rodek22"]["control"]


# ===========================================================================
print("=" * 78)
print("S7 -- EXACT TESTS  (cells are small; asymptotic p-values are optimistic)")
print("=" * 78)

ac, at = alleles(RODEK_C), alleles(RODEK_T)
tests = {
    "allelic (C vs A)":        [[ac["C"], ac["A"]], [at["C"], at["A"]]],
    "recessive (CC vs rest)":  [[RODEK_C["CC"], RODEK_C["AA"] + RODEK_C["AC"]],
                                [RODEK_T["CC"], RODEK_T["AA"] + RODEK_T["AC"]]],
    "dominant (CC+AC vs AA)":  [[RODEK_C["CC"] + RODEK_C["AC"], RODEK_C["AA"]],
                                [RODEK_T["CC"] + RODEK_T["AC"], RODEK_T["AA"]]],
}
print(f"{'contrast':<24} | {'chi2 p':>10} | {'Fisher p':>10} | {'mid-p':>10}")
for name, tbl in tests.items():
    t = np.array(tbl)
    _, p_chi, _, _ = stats.chi2_contingency(t, correction=False)
    odds, p_fisher = stats.fisher_exact(t)
    # mid-p correction: Fisher's exact is conservative
    tab = stats.hypergeom(t.sum(), t[:, 0].sum(), t[0].sum())
    k = t[0, 0]
    p_mid = 2 * min(tab.cdf(k) - 0.5 * tab.pmf(k),
                    tab.sf(k) + 0.5 * tab.pmf(k))
    print(f"{name:<24} | {p_chi:>10.2e} | {p_fisher:>10.2e} | {min(1,p_mid):>10.2e}")
print("\n  The allelic test treats 322 alleles from 161 people as independent")
print("  observations. The genotypic (recessive/dominant) tests do not, and")
print("  are the appropriate basis for inference.")


# ===========================================================================
print("\n" + "=" * 78)
print("S8 -- BAYESIAN REANALYSIS")
print("A p-value answers the wrong question. Under a prior calibrated to")
print("known common-variant effects in depression, what is the posterior?")
print("=" * 78)

y, se = or_2x2(ac["C"], ac["A"], at["C"], at["A"])
print(f"\n  observed: log OR = {y:.4f}  SE = {se:.4f}  (OR {exp(y):.3f})\n")
print(f"{'prior SD':>9} | {'prior 95% OR range':>22} | {'posterior OR':>13} | "
      f"{'shrinkage':>10} | {'BF10':>7}")
for tau in (0.05, 0.10, 0.20, 0.50):
    prec = 1 / tau ** 2 + 1 / se ** 2
    mu_post = (y / se ** 2) / prec
    se_post = sqrt(1 / prec)
    bf10 = (stats.norm.pdf(y, 0, sqrt(tau ** 2 + se ** 2))
            / stats.norm.pdf(y, 0, se))
    print(f"{tau:>9.2f} | {exp(-1.96*tau):>10.2f}-{exp(1.96*tau):<11.2f} | "
          f"{exp(mu_post):>13.3f} | {1 - mu_post/y:>9.0%} | {bf10:>7.2f}")

print("\n  A prior SD of 0.05 (95% of true effects between OR 0.91 and 1.10)")
print("  is the realistic one for a common variant in a polygenic disorder.")
print("  Under it the published OR of 1.57 shrinks to ~1.06, and the Bayes")
print("  factor is ~2 -- 'barely worth a bare mention' on Jeffreys' scale --")
print("  despite a p-value below 0.001.")


# ===========================================================================
print("\n" + "=" * 78)
print("S9 -- QUANTITATIVE BIAS ANALYSIS: COULD GENOTYPING ERROR EXPLAIN IT?")
print("Testing our own suggested explanation rather than asserting it.")
print("=" * 78)

n_case = sum(RODEK_C.values())
n_ctl = sum(RODEK_T.values())
true_dist = np.array([RODEK_T["AA"], RODEK_T["AC"], RODEK_T["CC"]]) / n_ctl
observed = np.array([RODEK_C["AA"], RODEK_C["AC"], RODEK_C["CC"]])


def called_counts(eps, n, true_p):
    """Per-allele mis-scoring: a true A allele is read as C with prob eps.

    Over-digestion (star activity) in PCR-RFLP shifts calls toward the cut
    (C) allele, which is the direction required to inflate CC in cases.
    """
    aa, ac, cc = true_p
    p_aa = aa * (1 - eps) ** 2
    p_ac = aa * 2 * eps * (1 - eps) + ac * (1 - eps)
    p_cc = aa * eps ** 2 + ac * eps + cc
    return np.array([p_aa, p_ac, p_cc]) * n


def gap(eps):
    return called_counts(eps, n_case, true_dist)[2] - observed[2]


eps_hat = optimize.brentq(gap, 1e-6, 0.95)
fitted = called_counts(eps_hat, n_case, true_dist)
print(f"\n  Null model: the true case distribution equals the control")
print(f"  distribution, and the apparent effect is mis-scoring alone.\n")
print(f"  Required per-allele mis-scoring rate: {eps_hat:.3f}  "
      f"({eps_hat*100:.1f}%)")
print(f"\n  {'genotype':<10} | {'observed':>9} | {'fitted':>9} | {'true (null)':>12}")
for i, g in enumerate(("AA", "AC", "CC")):
    print(f"  {g:<10} | {observed[i]:>9} | {fitted[i]:>9.1f} | "
          f"{true_dist[i]*n_case:>12.1f}")

chi_fit = float(((observed - fitted) ** 2 / fitted).sum())
print(f"\n  Goodness of fit, 2 df: chi2 = {chi_fit:.1f}, "
      f"p = {1 - stats.chi2.cdf(chi_fit, 2):.4f}")
print("  Even at this implausible rate the model does not reproduce the")
print("  observed distribution: it is fitted to CC and then under-predicts")
print("  AA by a third. No single mis-scoring rate generates the data.")

print(f"\n  Reported genotyping QC: 5% of samples re-genotyped, 100% concordance.")
print(f"  An error rate of {eps_hat*100:.0f}% would have produced roughly "
      f"{1-(1-eps_hat)**(2*8):.0%} discordance")
print(f"  across 8 re-genotyped samples, which was not observed.")
print("\n  CONCLUSION: differential genotyping error at any plausible rate")
print("  (1-5%) cannot generate the observed effect. The heterozygote deficit")
print("  reported in S3 is therefore NOT a sufficient explanation, and we")
print("  do not advance it as one. The anomaly localises to the case series,")
print("  but its source is more likely case selection, cryptic substructure,")
print("  or chance.")

# how much of the effect could a plausible error rate produce?
print(f"\n  {'error rate':>11} | {'apparent OR from a true null':>30}")
for eps in (0.01, 0.02, 0.05, 0.10):
    f = called_counts(eps, n_case, true_dist)
    fa = {"AA": f[0], "AC": f[1], "CC": f[2]}
    b, _ = or_2x2(fa["CC"], fa["AA"] + fa["AC"],
                  RODEK_T["CC"], RODEK_T["AA"] + RODEK_T["AC"])
    print(f"{eps:>11.0%} | {exp(b):>30.3f}")
print(f"  (observed recessive OR = 2.29)")


# ===========================================================================
print("\n" + "=" * 78)
print("S10 -- VIBRATION OF EFFECTS")
print("The range of odds ratios obtainable from one dataset by varying only")
print("defensible analytic choices.")
print("=" * 78)

# sex-stratified counts, Rodek Table 2 and Table 5
strata = {
    "all":   {"case": RODEK_C, "control": RODEK_T},
    "women": {"case": {"AA": 30, "AC": 51, "CC": 35},
              "control": {"AA": 73, "AC": 131, "CC": 34}},
    "men":   {"case": {"AA": 11, "AC": 20, "CC": 14},
              "control": {"AA": 81, "AC": 103, "CC": 40}},
}
contrasts = {
    "allelic":     lambda g: (alleles(g)["C"], alleles(g)["A"]),
    "recessive":   lambda g: (g["CC"], g["AA"] + g["AC"]),
    "dominant":    lambda g: (g["CC"] + g["AC"], g["AA"]),
    "overdominant": lambda g: (g["AC"], g["AA"] + g["CC"]),
    "CC vs AA":    lambda g: (g["CC"], g["AA"]),
}
rows = []
for (sname, s), (cname, fn) in product(strata.items(), contrasts.items()):
    a, b = fn(s["case"]); c, d = fn(s["control"])
    lo, se = or_2x2(a, b, c, d)
    p = 2 * (1 - stats.norm.cdf(abs(lo / se)))
    rows.append((sname, cname, exp(lo), p))

ors = [r[2] for r in rows]
sig = [r for r in rows if r[3] < 0.05]
print(f"\n  {len(rows)} defensible analyses of one dataset")
print(f"  odds ratios range from {min(ors):.2f} to {max(ors):.2f} "
      f"({max(ors)/min(ors):.1f}-fold)")
print(f"  {len(sig)}/{len(rows)} reach p < 0.05 unadjusted")
print(f"\n  {'stratum':<7} | {'contrast':<13} | {'OR':>6} | {'p':>9}")
for s, c, o, p in sorted(rows, key=lambda r: -r[2]):
    mark = " *" if p < 0.05 else ""
    print(f"  {s:<7} | {c:<13} | {o:>6.2f} | {p:>9.4f}{mark}")
print("\n  The reported result (recessive, all) is not the largest available.")
print("  The spread is what a reader needs in order to judge a single")
print("  selected estimate, and it is why the inheritance model must be")
print("  pre-specified rather than chosen by fit.")


# ===========================================================================
print("\n" + "=" * 78)
print("S11 -- SIMULATION: WHAT WOULD THIS LITERATURE LOOK LIKE UNDER A NULL?")
print("=" * 78)

rng = np.random.default_rng(20260207)
N_SIM = 40000

POWERS = [power_allelic(r["n_cases"], r["n_controls"], r["control_freq_C"],
                        1.0) for r in L]


def frac_significant_among_published(true_or, q_null, n_sim=N_SIM):
    """Conduct studies, suppress non-significant ones with prob (1-q_null),
    then look at the published record.

    Publication bias acts by making conducted studies invisible, so the
    quantity to model is the fraction significant AMONG PUBLISHED, not
    among conducted. Studies are drawn with replacement from the six
    located designs until six are published.
    """
    pw = np.array([power_allelic(r["n_cases"], r["n_controls"],
                                 r["control_freq_C"], true_or) for r in L])
    n_sig = np.zeros(n_sim, dtype=int)
    for _ in range(6):
        idx = rng.integers(0, len(L), n_sim)
        sig = rng.random(n_sim) < pw[idx]
        # resample until published: P(sig | published)
        p_pub_sig = pw[idx]
        p_pub_null = (1 - pw[idx]) * q_null
        p_sig_given_pub = p_pub_sig / (p_pub_sig + p_pub_null)
        n_sig += rng.random(n_sim) < p_sig_given_pub
    return n_sig


print("\n  Observed: 2 of 6 published studies significant (33%).")
print("  Publication bias is modelled as suppression of non-significant")
print("  findings, so the diagnostic quantity is the fraction significant")
print("  among PUBLISHED studies.\n")
print(f"  {'true OR':>8} | " + " | ".join(
    f"{f'q={q}':>9}" for q in (1.0, 0.5, 0.2, 0.1)))
print(f"  {'':>8} | " + " | ".join(f"{'':>9}" for _ in range(4)))
for true_or in (1.00, 1.02, 1.05, 1.10, 1.20):
    cells = []
    for q in (1.0, 0.5, 0.2, 0.1):
        ns = frac_significant_among_published(true_or, q)
        cells.append(f"{np.mean(ns >= 2):>8.1%}")
    print(f"  {true_or:>8.2f} | " + " | ".join(cells))
print("\n  q = probability a non-significant study is published.")
print("  Cells give P(>= 2 of 6 published studies significant).")
print("\n  Under a true null with complete reporting the observed pattern is")
print("  uncommon. Under a true null with 80-90% of negative findings")
print("  unpublished it is unremarkable. A literature of six small studies")
print("  cannot separate a small true effect from no effect plus selective")
print("  reporting; it can exclude a large one.")


print("S12 -- LEAVE-ONE-OUT AND CUMULATIVE META-ANALYSIS")
print("=" * 78)

mdd = ["rodek22", "thaweethee19"]


def allelic(sid):
    c, t = G[sid]["case"], G[sid]["control"]
    a, b = alleles(c), alleles(t)
    return or_2x2(a["C"], a["A"], b["C"], b["A"])


print("\n  Leave-one-out (k = 2, so each is a single study):")
for drop in mdd:
    keep = [s for s in mdd if s != drop]
    lo, se = allelic(keep[0])
    print(f"    omitting {G[drop]['meta']['study']:<18} -> "
          f"OR {exp(lo):.2f} ({exp(lo-1.96*se):.2f}-{exp(lo+1.96*se):.2f})")

print("\n  Cumulative by year:")
order = sorted(mdd, key=lambda s: G[s]["meta"]["year"])
ys, ses = [], []
for sid in order:
    lo, se = allelic(sid)
    ys.append(lo); ses.append(se)
    if len(ys) == 1:
        print(f"    {G[sid]['meta']['year']}  OR {exp(lo):.2f} "
              f"({exp(lo-1.96*se):.2f}-{exp(lo+1.96*se):.2f})   [k=1]")
    else:
        r = meta(ys, ses)
        print(f"    {G[sid]['meta']['year']}  OR {exp(r['mu_f']):.2f} "
              f"({exp(r['mu_f']-1.96*r['se_f']):.2f}-"
              f"{exp(r['mu_f']+1.96*r['se_f']):.2f})   [k={len(ys)}]")
print("\n  With k = 2 these are presented for completeness, not inference.")
print("  Neither can detect the drift a cumulative analysis is designed to")
print("  reveal; both are reported because their absence would be noted.")
