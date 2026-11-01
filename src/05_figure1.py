#!/usr/bin/env python3
"""Figure 1 for the manuscript: four-panel forensic summary."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from math import log, exp, sqrt

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
})
INK, ACC, MUT = "#1a1a1a", "#b2182b", "#8c8c8c"

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.32,
                      left=0.13, right=0.97, top=0.93, bottom=0.09)

# ---------------------------------------------------------------- A: forest
ax = fig.add_subplot(gs[0, 0])
rows = [
    ("Rodek 2022  MDD  PL*", 1.567, 1.215, 2.022, ACC),
    ("Thaweethee 2019  MDD  TH", 1.157, 0.751, 1.783, INK),
    ("Merk 2019  SZ  PL*", 1.080, 0.890, 1.310, INK),
    ("", None, None, None, None),
    ("Pooled, fixed effect", 1.449, 1.163, 1.805, MUT),
    ("Pooled, Hartung–Knapp", 1.419, 0.232, 8.670, ACC),
]
ys = list(range(len(rows)))[::-1]
for y, (lab, o, lo, hi, c) in zip(ys, rows):
    if o is None:
        continue
    ax.plot([lo, hi], [y, y], color=c, lw=1.1, solid_capstyle="butt")
    ax.plot([lo, lo], [y-.14, y+.14], color=c, lw=1.1)
    ax.plot([hi, hi], [y-.14, y+.14], color=c, lw=1.1)
    ax.plot(o, y, "s", color=c, ms=4.5, mec="white", mew=0.5)
ax.axvline(1, color=MUT, lw=0.7, ls=(0, (3, 3)))
ax.set_xscale("log")
ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8"])
ax.set_yticks(ys)
ax.set_yticklabels([r[0] for r in rows], fontsize=7)
ax.set_xlabel("Odds ratio, C allele (log scale)", fontsize=7.5)
ax.set_ylim(-0.6, len(rows) - 0.3)
ax.text(0.0, -0.20, "* shared control series (n = 462)",
        transform=ax.transAxes, fontsize=6, color=MUT, style="italic")
ax.set_title("a   Allelic association, rs4354668", loc="left",
             fontweight="bold", fontsize=8.5, pad=6)

# ---------------------------------------------------------------- B: het deficit
ax = fig.add_subplot(gs[0, 1])
labs = ["MDD PL\n161", "SZ PL\n393", "Ctrl PL\n462", "MDD TH\n100", "Ctrl TH\n100"]
F = [0.116, 0.026, -0.044, 0.022, -0.109]
lo = [-0.039, -0.075, -0.136, -0.178, -0.290]
hi = [0.273, 0.121, 0.048, 0.207, 0.080]
cols = [ACC] + [INK]*4
xs = np.arange(len(F))
for x, f, l, h, c in zip(xs, F, lo, hi, cols):
    ax.plot([x, x], [l, h], color=c, lw=1.1)
    ax.plot(x, f, "o", color=c, ms=5, mec="white", mew=0.6)
ax.axhline(0, color=MUT, lw=0.7, ls=(0, (3, 3)))
ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=6.5)
ax.set_xlim(-0.6, 4.6)
ax.set_ylabel("Wright's $F$ (heterozygote deficit)", fontsize=7.5)
ax.set_title("b   Genotyping quality control", loc="left",
             fontweight="bold", fontsize=8.5, pad=6)

# ---------------------------------------------------------------- C: excess sig
ax = fig.add_subplot(gs[1, 0])
ORs = [1.02, 1.05, 1.10, 1.20, 1.45]
E = [0.23, 0.41, 0.93, 2.19, 4.29]
ax.plot(ORs, E, "-o", color=INK, lw=1.2, ms=4, label="Expected (from power)")
ax.axhline(2, color=ACC, lw=1.4, label="Observed = 2")
ax.fill_between([1.02, 1.10], 0, 5, color=ACC, alpha=0.07)
ax.text(1.06, 3.4, "plausible\ntrue effect", fontsize=6, color=ACC,
        ha="center", va="top")
ax.annotate("$p$ = 0.019", xy=(1.02, 0.23), xytext=(1.05, 1.35), fontsize=6.5,
            color=ACC, arrowprops=dict(arrowstyle="-", color=ACC, lw=0.7))
ax.set_xlabel("Assumed true odds ratio", fontsize=7.5)
ax.set_ylabel("Number of significant findings", fontsize=7.5)
ax.set_ylim(0, 5)
ax.legend(frameon=False, fontsize=6.5, loc="lower right")
ax.set_title("c   Excess significance across 6 studies", loc="left",
             fontweight="bold", fontsize=8.5, pad=6)

# ---------------------------------------------------------------- D: AUC
ax = fig.add_subplot(gs[1, 1])
def auc(q, OR):
    pop = np.array([(1-q)**2, 2*q*(1-q), q**2])
    case = pop*np.array([1., OR, OR**2]); case /= case.sum()
    ctrl = pop.copy()
    a = 0.
    for i in range(3):
        for j in range(3):
            a += case[i]*ctrl[j]*(1. if i > j else 0.5 if i == j else 0.)
    return a
gridOR = np.linspace(1.0, 3.0, 200)
ax.plot(gridOR, [auc(0.41, o) for o in gridOR], color=INK, lw=1.4)
ax.axhline(0.75, color=ACC, lw=1.0, ls=(0, (4, 3)))
ax.text(1.05, 0.757, "clinical usability threshold", fontsize=6, color=ACC)
ax.axhline(0.5, color=MUT, lw=0.7, ls=(0, (3, 3)))
ax.text(2.72, 0.507, "chance", fontsize=6, color=MUT)
ax.plot(1.567, auc(0.41, 1.567), "o", color=ACC, ms=5, mec="white", mew=0.6)
ax.annotate("published OR 1.57\nAUC 0.584", xy=(1.567, auc(0.41, 1.567)),
            xytext=(1.95, 0.565), fontsize=6.5, color=ACC,
            arrowprops=dict(arrowstyle="->", color=ACC, lw=0.7))
ax.set_xlabel("Per-allele odds ratio", fontsize=7.5)
ax.set_ylabel("AUC (single variant)", fontsize=7.5)
ax.set_ylim(0.45, 0.82)
ax.set_title("d   Ceiling on diagnostic discrimination", loc="left",
             fontweight="bold", fontsize=8.5, pad=6)

fig.savefig("results/Figure1.png", dpi=400, bbox_inches="tight",
            facecolor="white")
fig.savefig("results/Figure1.pdf", bbox_inches="tight", facecolor="white")
print("Figure 1 written (PNG 400 dpi + vector PDF)")
