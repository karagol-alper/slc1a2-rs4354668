#!/usr/bin/env python3
"""Figure 3: PRISMA 2020 flow of studies.

Counts are enumerated from the named records in Supplementary Table S1 and
asserted here, so the diagram cannot drift from the study list.

    python src/09_figure3_prisma.py   -> results/Figure3.{png,pdf}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
INK, ACC, MUT = "#1a1a1a", "#b2182b", "#8c8c8c"

# ---------------------------------------------------------------- the records
# Every record is named, so the arithmetic is auditable.
ELIGIBLE = ["Rodek 2022", "Thaweethee 2019"]

EXCLUDED = {
    "Schizophrenia case-control": [
        "Merk 2019", "Zhang 2015", "Nagai 2009", "Deng 2004", "Wang 2022"],
    "Bipolar / schizophrenia, other designs": [
        "Dallaspezia 2012", "Poletti 2014a", "Poletti 2014b",
        "Spangaro 2012", "Spangaro 2014", "Fiorentino 2015"],
    "Suicidal behaviour": ["Murphy 2011", "Choi 2018"],
    "Quantitative endophenotype": ["Veldic 2019"],
    "Healthy volunteers, personality trait": ["Matsumoto 2007"],
    "Post-mortem expression, no case-control": ["Kerdsan 2012"],
    "Parkinson's disease": ["Appenzeller 2013"],
}
# excluded from the primary analysis but retained for the cross-disorder
# and excess-significance analyses
RETAINED = ["Merk 2019", "Zhang 2015", "Nagai 2009", "Deng 2004"]

n_exc = sum(len(v) for v in EXCLUDED.values())
n_elig = len(ELIGIBLE)
n_assessed = n_exc + n_elig
assert n_assessed == 19 and n_exc == 17 and n_elig == 2, "flow must sum"

fig, ax = plt.subplots(figsize=(7.2, 5.4))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def box(x, y, w, h, lines, bold_first=True, edge=INK, fill="white", fs=6.8):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.5,rounding_size=1.0",
        linewidth=0.9, edgecolor=edge, facecolor=fill, zorder=2))
    n = len(lines)
    for i, ln in enumerate(lines):
        ax.text(x + w / 2, y + h - (h / (n + 1)) * (i + 1), ln,
                ha="center", va="center", fontsize=fs, zorder=3,
                fontweight="bold" if (i == 0 and bold_first) else "normal",
                color=INK)


def arrow(x1, y1, x2, y2, col=INK):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), zorder=1,
                                 arrowstyle="-|>", mutation_scale=8,
                                 linewidth=0.9, color=col))


for yy, lab in ((85, "Identification"), (57, "Screening"), (22, "Included")):
    ax.text(2, yy, lab, rotation=90, ha="center", va="center",
            fontsize=7, color=MUT, fontweight="bold")

# --- identification
box(7, 80, 51, 13, [
    "Records identified",
    "8 structured search strategies, PubMed/MEDLINE and",
    "Google Scholar, plus hand-search of reference lists"])
arrow(32.5, 80, 32.5, 71)

# --- eligibility
box(7, 56, 51, 15, [
    f"Reports assessed for eligibility  (n = {n_assessed})",
    "Full text examined; every record named in Table S1",
    "",
    "Eligibility: case-control design, MDD diagnosis,",
    "genotype counts reported"], fs=6.6)
arrow(58, 63.5, 64, 63.5)
arrow(32.5, 56, 32.5, 36)

# --- exclusions
lines = [f"Excluded from primary analysis  (n = {n_exc})", ""]
for reason, recs in EXCLUDED.items():
    lines.append(f"{reason}  ({len(recs)})")
box(64, 46, 35, 27, lines, fs=6.0)

# --- included
box(7, 22, 51, 14, [
    f"Included in primary meta-analysis  (n = {n_elig})",
    "Rodek 2022 (Poland)   161 cases / 462 controls",
    "Thaweethee 2019 (Thailand)   100 / 100",
    "total  261 cases / 562 controls"], edge=ACC, fs=6.6)
arrow(58, 29, 64, 29, col=MUT)
box(64, 22, 35, 14, [
    "Retained for cross-disorder",
    "and bias analyses  (n = 4)",
    "Merk 2019, Zhang 2015,",
    "Nagai 2009, Deng 2004"],
    bold_first=False, edge=MUT, fill="#f7f7f7", fs=6.2)

ax.text(7, 15,
        "No eligible study of $\\it{SLC1A3}$ and MDD diagnosis was identified.",
        fontsize=7, color=ACC, ha="left")
ax.text(7, 10.5,
        "Title and abstract screening was iterative rather than a single indexed "
        "pass, so counts are reported from full-text\nassessment onward; all "
        "assessed records are named, and the search strategies are listed in "
        "Table S1.",
        fontsize=6.0, color=MUT, ha="left", va="top", linespacing=1.6)

out = Path(__file__).resolve().parents[1] / "results"
out.mkdir(exist_ok=True)
fig.savefig(out / "Figure3.png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(out / "Figure3.pdf", bbox_inches="tight", facecolor="white")
print(f"Figure 3 written to {out}  (assessed {n_assessed} = "
      f"{n_exc} excluded + {n_elig} included)")
