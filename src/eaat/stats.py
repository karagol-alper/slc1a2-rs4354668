"""Statistical primitives for the rs4354668 re-appraisal.

Every function here is used by at least two analysis scripts. Nothing in this
module reads data; nothing prints. See eaat.data for loading.
"""
from __future__ import annotations

from math import exp, log, sqrt

import numpy as np
from scipy import stats

__all__ = [
    "alleles", "hwe", "wright_f", "or_2x2", "fmt_or",
    "meta", "power_allelic", "n_required", "fprp",
    "fst_two_pop", "auc_additive", "tost",
]


# --------------------------------------------------------------------------
# genotype helpers
# --------------------------------------------------------------------------
def alleles(g: dict) -> dict:
    """Allele counts from a genotype dict {'AA','AC','CC'}."""
    return {"A": 2 * g["AA"] + g["AC"], "C": 2 * g["CC"] + g["AC"]}


def hwe(g: dict) -> tuple[float, float]:
    """Pearson chi-square goodness-of-fit to Hardy-Weinberg (1 df)."""
    n = sum(g.values())
    p = alleles(g)["A"] / (2 * n)
    expected = np.array([p ** 2, 2 * p * (1 - p), (1 - p) ** 2]) * n
    observed = np.array([g["AA"], g["AC"], g["CC"]])
    chi = float(((observed - expected) ** 2 / expected).sum())
    return chi, float(1 - stats.chi2.cdf(chi, 1))


def wright_f(g: dict, n_boot: int = 4000, seed: int = 7):
    """Wright's F (heterozygote deficit) with a multinomial bootstrap CI.

    F > 0 means fewer heterozygotes than Hardy-Weinberg predicts: the
    signature of under-digestion in PCR-RFLP, which miscalls A/C as a
    homozygote. HWE testing alone is insensitive to this.
    """
    n = sum(g.values())
    p = alleles(g)["A"] / (2 * n)
    exp_het = 2 * p * (1 - p) * n
    f = 1 - g["AC"] / exp_het

    counts = np.array([g["AA"], g["AC"], g["CC"]])
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n, counts / n, size=n_boot)
    fs = []
    for d in draws:
        pp = (2 * d[0] + d[1]) / (2 * n)
        eh = 2 * pp * (1 - pp) * n
        if eh > 0:
            fs.append(1 - d[1] / eh)
    lo, hi = np.percentile(fs, [2.5, 97.5])
    return float(f), (float(lo), float(hi))


# --------------------------------------------------------------------------
# effect estimates
# --------------------------------------------------------------------------
def or_2x2(a, b, c, d) -> tuple[float, float]:
    """log OR and Woolf SE. a/b = exposed/unexposed cases, c/d = controls.

    Applies the Haldane-Anscombe 0.5 correction only if a cell is empty.
    """
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return log((a / b) / (c / d)), sqrt(1 / a + 1 / b + 1 / c + 1 / d)


def fmt_or(lor: float, se: float, dp: int = 2) -> str:
    return (f"{exp(lor):.{dp}f} "
            f"({exp(lor - 1.96 * se):.{dp}f}-{exp(lor + 1.96 * se):.{dp}f})")


def meta(y, s) -> dict:
    """Fixed (IV) and random (DerSimonian-Laird) pooling, plus Hartung-Knapp.

    At small k the DerSimonian-Laird variance is severely anticonservative,
    so the Hartung-Knapp t-interval is the result to report.
    """
    y, s = np.asarray(y, float), np.asarray(s, float)
    w = 1 / s ** 2
    mu_f = float((w * y).sum() / w.sum())
    se_f = float(sqrt(1 / w.sum()))

    k = len(y)
    df = k - 1
    Q = float((w * (y - mu_f) ** 2).sum())
    pQ = float(1 - stats.chi2.cdf(Q, df)) if df > 0 else float("nan")
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0

    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / C)
    wr = 1 / (s ** 2 + tau2)
    mu_r = float((wr * y).sum() / wr.sum())
    se_r = float(sqrt(1 / wr.sum()))

    if df > 0:
        se_hk = float(sqrt(((wr * (y - mu_r) ** 2).sum() / df) / wr.sum()))
        t_crit = float(stats.t.ppf(0.975, df))
    else:
        se_hk = t_crit = float("nan")

    return dict(k=k, mu_f=mu_f, se_f=se_f, mu_r=mu_r, se_r=se_r,
                se_hk=se_hk, t=t_crit, Q=Q, pQ=pQ, I2=I2, tau2=tau2)


# --------------------------------------------------------------------------
# power, credibility, bias
# --------------------------------------------------------------------------
def power_allelic(n_case, n_ctrl, maf_ctrl, OR, alpha=0.05) -> float:
    """Power of the allelic chi-square test (2N alleles per group)."""
    p0 = maf_ctrl
    p1 = OR * p0 / (1 - p0 + OR * p0)
    n1, n0 = 2 * n_case, 2 * n_ctrl
    pbar = (n1 * p1 + n0 * p0) / (n1 + n0)
    se0 = sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n0))
    se1 = sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    z = stats.norm.ppf(1 - alpha / 2)
    return float(1 - stats.norm.cdf((z * se0 - abs(p1 - p0)) / se1))


def n_required(maf, OR, alpha, power=0.80, ratio=1.0) -> int:
    """Cases needed for `power`, by bisection."""
    lo, hi = 10, 50_000_000
    while hi - lo > 10:
        mid = (lo + hi) // 2
        if power_allelic(mid, int(mid * ratio), maf, OR, alpha) < power:
            lo = mid
        else:
            hi = mid
    return hi


def fprp(p_obs, power, prior) -> float:
    """False-positive report probability (Wacholder et al. 2004).

    `power` must be evaluated at the PLAUSIBLE true effect. Using the
    observed effect is the winner's-curse error this quantity exists to
    expose.
    """
    return (p_obs * (1 - prior)) / (p_obs * (1 - prior) + power * prior)


def fst_two_pop(p1: float, p2: float) -> float:
    """Wright's FST between two populations at one locus."""
    pbar = (p1 + p2) / 2
    return float(((p1 - p2) ** 2 / 4) / (pbar * (1 - pbar)))


# --------------------------------------------------------------------------
# clinical utility
# --------------------------------------------------------------------------
def auc_additive(risk_freq: float, OR: float) -> float:
    """AUC of a single biallelic variant under a log-additive model."""
    q = risk_freq
    pop = np.array([(1 - q) ** 2, 2 * q * (1 - q), q ** 2])
    case = pop * np.array([1.0, OR, OR ** 2])
    case /= case.sum()
    ctrl = pop / pop.sum()
    a = 0.0
    for i in range(3):
        for j in range(3):
            a += case[i] * ctrl[j] * (1.0 if i > j else 0.5 if i == j else 0.0)
    return float(a)


# --------------------------------------------------------------------------
# equivalence
# --------------------------------------------------------------------------
def tost(beta: float, se: float, bound_or: float = 1.02, alpha: float = 0.05):
    """Two one-sided tests for equivalence against an odds-ratio bound.

    Returns (p_TOST, equivalent, tightest_bound_supported). A significant
    result is a POSITIVE finding of absence, not a failure to reject.
    """
    d = log(bound_or)
    p_lo = float(1 - stats.norm.cdf((beta + d) / se))
    p_hi = float(stats.norm.cdf((beta - d) / se))
    p = max(p_lo, p_hi)
    z = stats.norm.ppf(1 - alpha)
    ci = (beta - z * se, beta + z * se)
    equivalent = ci[0] > -d and ci[1] < d
    tightest = exp(max(abs(ci[0]), abs(ci[1])))
    return p, equivalent, tightest
