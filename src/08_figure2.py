#!/usr/bin/env python3
"""Figure 2: the extended analyses (Bayesian, bias analysis, vibration,
simulation). Values are computed here, not copied from the tables.

    python src/08_figure2.py      -> results/Figure2.{png,pdf}
"""
import sys
from itertools import product
from math import exp, sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eaat.data import load_genotypes, load_ledger
from eaat.stats import alleles, or_2x2, power_allelic

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
})
INK, ACC, MUT = "#1a1a1a", "#b2182b", "#8c8c8c"

G, L = load_genotypes(), load_ledger()
CS, CT = G["rodek22"]["case"], G["rodek22"]["control"]
ac, at = alleles(CS), alleles(CT)
y, se = or_2x2(ac["C"], ac["A"], at["C"], at["A"])

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.34,
                      left=0.11, right=0.97, top=0.93, bottom=0.085)

# ------------------------------------------------------- a: Bayesian shrinkage
ax = fig.add_subplot(gs[0, 0])
taus = np.array([0.05, 0.10, 0.20, 0.50])
post = [(y / se**2) / (1 / t**2 + 1 / se**2) for t in taus]
bf = [stats.norm.pdf(y, 0, sqrt(t**2 + se**2)) / stats.norm.pdf(y, 0, se)
      for t in taus]
ypos = np.arange(len(taus))[::-1]
ax.axvline(exp(y), color=ACC, lw=1.2, ls=(0, (4, 3)))
ax.text(exp(y), len(taus) - 0.35, " published\n OR 1.57", fontsize=6.3,
        color=ACC, va="top", ha="left")
ax.axvline(1, color=MUT, lw=0.7, ls=(0, (3, 3)))
for yy, t, pm, b in zip(ypos, taus, post, bf):
    ax.plot([1, exp(pm)], [yy, yy], color=MUT, lw=0.8, zorder=1)
    ax.plot(exp(pm), yy, "o", color=INK if t > 0.05 else ACC,
            ms=6 if t == 0.05 else 4.5, mec="white", mew=0.6, zorder=3)
    ax.text(exp(pm), yy + 0.26, f"BF {b:.1f}", fontsize=6, ha="center",
            color=ACC if t == 0.05 else MUT)
ax.set_yticks(ypos)
ax.set_yticklabels([f"{t:.2f}" for t in taus], fontsize=7)
ax.set_ylabel("Prior SD on log OR", fontsize=7.5)
ax.set_xlabel("Posterior odds ratio", fontsize=7.5)
ax.set_xlim(0.95, 1.75)
ax.set_ylim(-0.6, len(taus) - 0.1)
ax.set_title("a   Posterior shrinkage of the published estimate",
             loc="left", fontweight="bold", fontsize=8.5, pad=6)

# ------------------------------------------- b: quantitative bias analysis
ax = fig.add_subplot(gs[0, 1])
n_case = sum(CS.values())
true_p = np.array([CT["AA"], CT["AC"], CT["CC"]]) / sum(CT.values())


def called(eps):
    aa, acf, cc = true_p
    return np.array([aa * (1 - eps) ** 2,
                     aa * 2 * eps * (1 - eps) + acf * (1 - eps),
                     aa * eps ** 2 + acf * eps + cc]) * n_case


grid = np.linspace(0, 0.30, 200)
ors = []
for e in grid:
    f = called(e)
    b, _ = or_2x2(max(f[2], 1e-9), f[0] + f[1], CT["CC"], CT["AA"] + CT["AC"])
    ors.append(exp(b))
obs_or = exp(or_2x2(CS["CC"], CS["AA"] + CS["AC"],
                    CT["CC"], CT["AA"] + CT["AC"])[0])
eps_hat = optimize.brentq(lambda e: called(e)[2] - CS["CC"], 1e-6, 0.95)

ax.plot(grid * 100, ors, color=INK, lw=1.4)
ax.axhline(obs_or, color=ACC, lw=1.2, ls=(0, (4, 3)))
ax.text(0.6, obs_or + 0.06, f"observed OR {obs_or:.2f}", fontsize=6.3, color=ACC)
ax.axvspan(1, 5, color=MUT, alpha=0.13)
ax.text(3, 0.75, "plausible\nerror rate", fontsize=6, color=MUT,
        ha="center", va="bottom")
ax.plot(eps_hat * 100, obs_or, "o", color=ACC, ms=5, mec="white", mew=0.6)
ax.annotate(f"{eps_hat*100:.1f}% required", xy=(eps_hat * 100, obs_or),
            xytext=(14, 1.55), fontsize=6.3, color=ACC,
            arrowprops=dict(arrowstyle="->", color=ACC, lw=0.7))
ax.set_xlabel("Per-allele mis-scoring rate (%)", fontsize=7.5)
ax.set_ylabel("Apparent recessive OR from a true null", fontsize=7.5)
ax.set_ylim(0.7, 2.6)
ax.set_xlim(0, 30)
ax.set_title("b   Genotyping error cannot explain the finding",
             loc="left", fontweight="bold", fontsize=8.5, pad=6)

# ------------------------------------------------- c: vibration of effects
ax = fig.add_subplot(gs[1, 0])
strata = {
    "all": (CS, CT),
    "women": ({"AA": 30, "AC": 51, "CC": 35}, {"AA": 73, "AC": 131, "CC": 34}),
    "men": ({"AA": 11, "AC": 20, "CC": 14}, {"AA": 81, "AC": 103, "CC": 40}),
}
contrasts = {
    "allelic": lambda g: (alleles(g)["C"], alleles(g)["A"]),
    "recessive": lambda g: (g["CC"], g["AA"] + g["AC"]),
    "dominant": lambda g: (g["CC"] + g["AC"], g["AA"]),
    "overdom.": lambda g: (g["AC"], g["AA"] + g["CC"]),
    "CC vs AA": lambda g: (g["CC"], g["AA"]),
}
rows = []
for (sn, (c_, t_)), (cn, fn) in product(strata.items(), contrasts.items()):
    a, b = fn(c_); c, d = fn(t_)
    lo, s_ = or_2x2(a, b, c, d)
    rows.append((exp(lo), 2 * (1 - stats.norm.cdf(abs(lo / s_))),
                 f"{sn}/{cn}", sn == "all" and cn == "recessive"))
rows.sort(key=lambda r: r[0])
xs = np.arange(len(rows))
for i, (o, p, lab, reported) in enumerate(rows):
    col = ACC if reported else (INK if p < 0.05 else MUT)
    # open circles for non-significant: edge must be the colour, not white,
    # or the marker disappears against the page
    ax.plot(i, o, "o", ms=7 if reported else 4.5, zorder=3,
            mfc=col if p < 0.05 else "white",
            mec=col, mew=0.9 if p >= 0.05 else 0.6)
ax.axhline(1, color=MUT, lw=0.7, ls=(0, (3, 3)))
rep_i = [i for i, r in enumerate(rows) if r[3]][0]
ax.annotate("reported\nresult", xy=(rep_i, rows[rep_i][0]),
            xytext=(rep_i - 4.2, 2.95), fontsize=6.3, color=ACC, ha="center",
            arrowprops=dict(arrowstyle="->", color=ACC, lw=0.7))
ax.set_xticks(xs)
ax.set_xticklabels([r[2] for r in rows], rotation=90, fontsize=5.3)
ax.set_ylabel("Odds ratio", fontsize=7.5)
ax.set_yscale("log")
ax.set_yticks([0.5, 1, 2, 3])
ax.set_yticklabels(["0.5", "1", "2", "3"])
ax.minorticks_off()
ax.set_ylim(0.45, 3.6)
ax.set_title("c   Vibration of effects across 15 analyses",
             loc="left", fontweight="bold", fontsize=8.5, pad=6)
ax.text(0.02, 0.96, "filled = $p$ < 0.05", transform=ax.transAxes,
        fontsize=6, color=MUT, va="top")

# ------------------------------------------------------------ d: simulation
ax = fig.add_subplot(gs[1, 1])
rng = np.random.default_rng(20260207)
N = 20000
true_ors = np.array([1.00, 1.02, 1.05, 1.10, 1.20])
styles = {1.0: (INK, "-"), 0.5: (INK, (0, (5, 2))),
          0.2: (ACC, (0, (3, 2))), 0.1: (ACC, "-")}
for q, (col, ls) in styles.items():
    vals = []
    for t_or in true_ors:
        pw = np.array([power_allelic(r["n_cases"], r["n_controls"],
                                     r["control_freq_C"], t_or) for r in L])
        n_sig = np.zeros(N, dtype=int)
        for _ in range(6):
            i = rng.integers(0, len(L), N)
            n_sig += rng.random(N) < pw[i] / (pw[i] + (1 - pw[i]) * q)
        vals.append(np.mean(n_sig >= 2) * 100)
    ax.plot(true_ors, vals, color=col, ls=ls, lw=1.3, marker="o", ms=3.2,
            label=f"$q$ = {q}")
ax.set_xlabel("True odds ratio", fontsize=7.5)
ax.set_ylabel("P(≥ 2 of 6 significant)  %", fontsize=7.5)
ax.legend(frameon=False, fontsize=6.2, loc="upper left",
          title="published if null", title_fontsize=6.2)
ax.set_ylim(-3, 103)
ax.set_title("d   Observed record under selective reporting",
             loc="left", fontweight="bold", fontsize=8.5, pad=6)

out = Path(__file__).resolve().parents[1] / "results"
out.mkdir(exist_ok=True)
fig.savefig(out / "Figure2.png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(out / "Figure2.pdf", bbox_inches="tight", facecolor="white")
print(f"Figure 2 written to {out}  (400 dpi PNG + vector PDF)")
