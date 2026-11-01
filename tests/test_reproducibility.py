"""Every number published in the manuscript and supplementary, asserted.

If a refactor, a dependency bump or an edit to data/genotype_counts.csv
changes any reported value, this suite fails. That is the point: the paper
documents a literature in which numbers drifted from their sources, and it
should not repeat the error.

Run from the repository root:  pytest -q
"""
import sys
from math import log, sqrt
from pathlib import Path

import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from eaat.data import SHARED_CONTROLS, load_genotypes, load_ledger
from eaat.stats import (alleles, auc_additive, fprp, fst_two_pop, hwe, meta,
                        n_required, or_2x2, power_allelic, tost, wright_f)

G = load_genotypes()
L = load_ledger()
A = pytest.approx


# ---------------------------------------------------------------- data
def test_genotype_counts_match_published_tables():
    assert G["rodek22"]["case"] == {"AA": 41, "AC": 71, "CC": 49}
    assert G["rodek22"]["control"] == {"AA": 154, "AC": 234, "CC": 74}
    assert G["merk19"]["case"] == {"AA": 129, "AC": 188, "CC": 76}
    assert G["thaweethee19"]["case"] == {"AA": 8, "AC": 39, "CC": 53}
    assert G["thaweethee19"]["control"] == {"AA": 7, "AC": 47, "CC": 46}


def test_rodek_and_merk_share_one_control_series():
    """The finding the manuscript turns on. Must not silently diverge."""
    assert G["rodek22"]["control"] == G["merk19"]["control"]
    assert SHARED_CONTROLS == {"rodek22", "merk19"}


def test_allele_counts_match_source_papers():
    a = alleles(G["rodek22"]["case"])
    assert (a["A"], a["C"]) == (153, 169)
    a = alleles(G["rodek22"]["control"])
    assert (a["A"], a["C"]) == (542, 382)


# ------------------------------------------------------- per-study effects
def _allelic(sid):
    cs, ct = G[sid]["case"], G[sid]["control"]
    ac, at = alleles(cs), alleles(ct)
    return or_2x2(ac["C"], ac["A"], at["C"], at["A"])


@pytest.mark.parametrize("sid,odds,lo,hi", [
    ("rodek22", 1.57, 1.21, 2.02),
    ("thaweethee19", 1.16, 0.75, 1.78),
    ("merk19", 1.08, 0.89, 1.31),
])
def test_allelic_odds_ratios(sid, odds, lo, hi):
    from math import exp
    b, se = _allelic(sid)
    assert exp(b) == A(odds, abs=0.005)
    assert exp(b - 1.96 * se) == A(lo, abs=0.005)
    assert exp(b + 1.96 * se) == A(hi, abs=0.005)


def test_rodek_recessive_odds_ratio():
    from math import exp
    cs, ct = G["rodek22"]["case"], G["rodek22"]["control"]
    b, se = or_2x2(cs["CC"], cs["AA"] + cs["AC"], ct["CC"], ct["AA"] + ct["AC"])
    assert exp(b) == A(2.29, abs=0.01)


def test_merk_schizophrenia_is_null():
    """Same lab, same controls, 2.4x more patients."""
    cs, ct = G["merk19"]["case"], G["merk19"]["control"]
    b, se = or_2x2(cs["CC"], cs["AA"] + cs["AC"], ct["CC"], ct["AA"] + ct["AC"])
    p = 2 * (1 - stats.norm.cdf(abs(b / se)))
    assert p > 0.05 and p == A(0.20, abs=0.01)


# --------------------------------------------------------------- pooling
def test_pooled_allelic_fixed_and_hartung_knapp():
    from math import exp
    ys, ses = zip(*(_allelic(s) for s in ("rodek22", "thaweethee19")))
    r = meta(ys, ses)
    assert exp(r["mu_f"]) == A(1.45, abs=0.005)
    assert exp(r["mu_f"] - 1.96 * r["se_f"]) == A(1.16, abs=0.005)
    assert exp(r["mu_f"] + 1.96 * r["se_f"]) == A(1.81, abs=0.005)
    lo = exp(r["mu_r"] - r["t"] * r["se_hk"])
    hi = exp(r["mu_r"] + r["t"] * r["se_hk"])
    assert (lo, hi) == (A(0.23, abs=0.01), A(8.67, abs=0.05))
    assert r["Q"] == A(1.40, abs=0.01)
    assert r["I2"] == A(28.8, abs=0.5)


def test_hartung_knapp_is_wider_than_dersimonian_laird():
    """The methodological claim: DL is anticonservative at k=2."""
    ys, ses = zip(*(_allelic(s) for s in ("rodek22", "thaweethee19")))
    r = meta(ys, ses)
    assert r["t"] * r["se_hk"] > 1.96 * r["se_r"]


# ----------------------------------------------------------- forensic checks
def test_heterozygote_deficit_largest_in_the_positive_series():
    f_mdd, (lo, hi) = wright_f(G["rodek22"]["case"])
    f_sz, _ = wright_f(G["merk19"]["case"])
    f_ctl, _ = wright_f(G["rodek22"]["control"])
    assert f_mdd == A(0.116, abs=0.001)
    assert f_sz == A(0.026, abs=0.001)
    assert f_ctl == A(-0.044, abs=0.001)
    assert f_mdd > f_sz > f_ctl
    assert lo < 0 < hi, "CI must include zero; the claim is not probative"


def test_all_series_pass_hardy_weinberg():
    for sid, rec in G.items():
        for grp in ("case", "control"):
            assert hwe(rec[grp])[1] > 0.05, f"{sid}/{grp}"


def test_case_case_comparison():
    from math import exp
    mdd, sz = G["rodek22"]["case"], G["merk19"]["case"]
    b, se = or_2x2(mdd["CC"], mdd["AA"] + mdd["AC"],
                   sz["CC"], sz["AA"] + sz["AC"])
    p = 2 * (1 - stats.norm.cdf(abs(b / se)))
    assert exp(b) == A(1.82, abs=0.01)
    assert p == A(0.005, abs=0.001)


def test_no_sex_by_genotype_interaction():
    from math import exp
    lf, sf = or_2x2(35, 81, 34, 204)
    lm, sm = or_2x2(14, 31, 40, 184)
    d, sd = lf - lm, sqrt(sf ** 2 + sm ** 2)
    assert exp(lf) == A(2.59, abs=0.01)
    assert exp(lm) == A(2.08, abs=0.01)
    assert exp(d) == A(1.25, abs=0.01)
    assert 2 * (1 - stats.norm.cdf(abs(d / sd))) == A(0.63, abs=0.01)


def test_ancestry_divergence():
    p1 = alleles(G["rodek22"]["control"])["C"] / 924
    p2 = alleles(G["thaweethee19"]["control"])["C"] / 200
    assert fst_two_pop(p1, p2) == A(0.080, abs=0.001)


# ------------------------------------------------------------ power, bias
@pytest.mark.parametrize("sid,OR,expected", [
    ("rodek22", 1.05, 5.7), ("rodek22", 1.02, 3.5), ("rodek22", 1.57, 93.5),
    ("thaweethee19", 1.05, 4.1), ("merk19", 1.05, 7.2),
])
def test_power(sid, OR, expected):
    rec = G[sid]
    maf = alleles(rec["control"])["C"] / (2 * sum(rec["control"].values()))
    pw = power_allelic(sum(rec["case"].values()),
                       sum(rec["control"].values()), maf, OR) * 100
    assert pw == A(expected, abs=0.1)


def test_required_sample_sizes():
    assert n_required(0.41, 1.05, 0.05) == A(13574, rel=0.01)
    assert n_required(0.41, 1.05, 5e-8) == A(68488, rel=0.01)


def test_fprp():
    maf = alleles(G["rodek22"]["control"])["C"] / 924
    pw = power_allelic(161, 462, maf, 1.05, alpha=0.001)
    assert fprp(0.001, pw, 0.05) == A(0.914, abs=0.005)
    pw = power_allelic(161, 462, maf, 1.20, alpha=0.001)
    assert fprp(0.001, pw, 0.05) == A(0.388, abs=0.005)


@pytest.mark.parametrize("OR,E,p", [
    (1.02, 0.23, 0.019), (1.05, 0.41, 0.059),
    (1.10, 0.93, 0.237), (1.20, 2.19, 0.707), (1.45, 4.29, 0.991),
])
def test_excess_significance(OR, E, p):
    O = sum(1 for r in L if r["significant"])
    e = sum(power_allelic(r["n_cases"], r["n_controls"],
                          r["control_freq_C"], OR) for r in L)
    assert e == A(E, abs=0.01)
    assert stats.binomtest(O, len(L), e / len(L),
                           alternative="greater").pvalue == A(p, abs=0.001)


def test_ledger_totals():
    assert sum(r["n_cases"] for r in L) == 2071
    assert sum(r["n_cases"] for r in L if r["phenotype"] == "MDD") == 261
    assert sum(r["n_controls"] for r in L if r["phenotype"] == "MDD") == 562


def test_sign_test_cannot_reach_significance_at_k4():
    dirs = [r["direction"] for r in L if r["direction"] != 0]
    p = stats.binomtest(sum(1 for d in dirs if d > 0), len(dirs), 0.5).pvalue
    assert p == A(0.125, abs=0.001) and p > 0.05


# -------------------------------------------------------- clinical utility
@pytest.mark.parametrize("OR,auc", [
    (1.05, 0.509), (1.10, 0.518), (1.20, 0.534),
    (1.57, 0.584), (2.00, 0.627), (3.00, 0.695),
])
def test_auc_ceiling(OR, auc):
    assert auc_additive(0.41, OR) == A(auc, abs=0.001)


def test_auc_at_published_effect_is_below_clinical_usability():
    assert auc_additive(0.41, 1.57) < 0.75


# ------------------------------------------------------------- equivalence
def test_tost_declares_equivalence_on_a_null_result():
    p, equiv, tight = tost(0.0012, 0.0009, bound_or=1.02)
    assert equiv and p < 1e-50 and tight == A(1.003, abs=0.001)


def test_tost_rejects_equivalence_at_the_published_effect():
    """Must not rubber-stamp. Fed the candidate-gene effect, it says no."""
    _, equiv, _ = tost(log(1.57), 0.0009, bound_or=1.02)
    assert not equiv


# ============================================================= extended (S7-S12)
def test_exact_tests_agree_with_asymptotic_in_direction():
    """Cells are small; Fisher must not overturn the asymptotic result."""
    from scipy import stats as st
    cs, ct = G["rodek22"]["case"], G["rodek22"]["control"]
    tbl = [[cs["CC"], cs["AA"] + cs["AC"]], [ct["CC"], ct["AA"] + ct["AC"]]]
    _, p_fisher = st.fisher_exact(tbl)
    assert p_fisher == A(1.35e-4, rel=0.05)
    assert p_fisher > st.chi2_contingency(tbl, correction=False)[1], \
        "Fisher's exact should be the more conservative of the two"


def test_bayesian_shrinkage_under_realistic_prior():
    """Published OR 1.57 collapses toward the null under a calibrated prior."""
    from math import exp
    cs, ct = G["rodek22"]["case"], G["rodek22"]["control"]
    a, b = alleles(cs), alleles(ct)
    y, se = or_2x2(a["C"], a["A"], b["C"], b["A"])
    tau = 0.05
    mu_post = (y / se ** 2) / (1 / tau ** 2 + 1 / se ** 2)
    assert exp(mu_post) == A(1.06, abs=0.01)
    bf10 = (stats.norm.pdf(y, 0, sqrt(tau ** 2 + se ** 2))
            / stats.norm.pdf(y, 0, se))
    assert bf10 == A(2.01, abs=0.05), "BF ~2 despite p < 0.001"
    assert bf10 < 3, "must remain below the 'substantial evidence' threshold"


def test_genotyping_error_cannot_explain_the_finding():
    """The paper's own suggested explanation, tested and found insufficient.

    Guards against the hypothesis hardening into a claim.
    """
    from scipy import optimize
    import numpy as np
    cs, ct = G["rodek22"]["case"], G["rodek22"]["control"]
    n_case = sum(cs.values())
    true_p = np.array([ct["AA"], ct["AC"], ct["CC"]]) / sum(ct.values())

    def called_cc(eps):
        aa, ac, cc = true_p
        return (aa * eps ** 2 + ac * eps + cc) * n_case

    eps = optimize.brentq(lambda e: called_cc(e) - cs["CC"], 1e-6, 0.95)
    assert eps == A(0.245, abs=0.005)
    assert eps > 0.15, "required error rate must stay implausibly high"


def test_plausible_error_rates_produce_negligible_bias():
    from math import exp
    import numpy as np
    cs, ct = G["rodek22"]["case"], G["rodek22"]["control"]
    n_case = sum(cs.values())
    aa, ac, cc = np.array([ct["AA"], ct["AC"], ct["CC"]]) / sum(ct.values())
    for eps, max_or in ((0.01, 1.05), (0.05, 1.25)):
        f_cc = (aa * eps ** 2 + ac * eps + cc) * n_case
        f_rest = n_case - f_cc
        b, _ = or_2x2(f_cc, f_rest, ct["CC"], ct["AA"] + ct["AC"])
        assert exp(b) < max_or


def test_vibration_of_effects_spans_a_fourfold_range():
    from math import exp
    from itertools import product
    strata = {
        "all": (G["rodek22"]["case"], G["rodek22"]["control"]),
        "women": ({"AA": 30, "AC": 51, "CC": 35},
                  {"AA": 73, "AC": 131, "CC": 34}),
        "men": ({"AA": 11, "AC": 20, "CC": 14},
                {"AA": 81, "AC": 103, "CC": 40}),
    }
    contrasts = [
        lambda g: (alleles(g)["C"], alleles(g)["A"]),
        lambda g: (g["CC"], g["AA"] + g["AC"]),
        lambda g: (g["CC"] + g["AC"], g["AA"]),
        lambda g: (g["AC"], g["AA"] + g["CC"]),
        lambda g: (g["CC"], g["AA"]),
    ]
    ors = []
    for (cs, ct), fn in product(strata.values(), contrasts):
        a, b = fn(cs); c, d = fn(ct)
        ors.append(exp(or_2x2(a, b, c, d)[0]))
    assert len(ors) == 15
    assert min(ors) == A(0.64, abs=0.01)
    assert max(ors) == A(2.59, abs=0.01)
    assert max(ors) / min(ors) > 3.5


def test_cumulative_meta_matches_pooled():
    """Cumulative analysis at k=2 must land on the pooled estimate."""
    from math import exp
    ys, ses = [], []
    for sid in ("thaweethee19", "rodek22"):
        cs, ct = G[sid]["case"], G[sid]["control"]
        a, b = alleles(cs), alleles(ct)
        y, se = or_2x2(a["C"], a["A"], b["C"], b["A"])
        ys.append(y); ses.append(se)
    assert exp(meta(ys, ses)["mu_f"]) == A(1.45, abs=0.005)
