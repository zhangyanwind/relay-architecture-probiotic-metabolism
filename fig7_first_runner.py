#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 7 — Consequences of being the first runner (Plan C, real computation)
=====================================================================
Replaces the old fig12_probiotic_intervention, which drew a bar chart of nine
hard-coded RCT sample sizes (GERD 120, BLa80 106, T2D 75, ...). Those values
were transcribed from a literature reading list and are not an output of this
study; the panel's claim ("interventions alter relay intermediates") could not
be tested with the pangenome data, and a bar chart of hand-entered numbers
invites rejection on sight. Under Plan C it is replaced by what the genomes can
actually establish: the quantitative consequence of a first-runner-only
metabolism.

(a) Pathway completion: who makes both halves of a relay, and who makes one.
(b) Orphan carriers - downstream enzymes without their own upstream supply.
(c) Minimum combination that completes the profiled pathways.
(d) Per-unit relay load.
"""
import sys, os, json
sys.path.insert(0, r"D:\Probiotics_Pangenome_v5\validation")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrowPatch
from paper_palette import (apply_style, PALETTE, UP_C, DOWN_C, BOND_C,
                           NEUTRAL, ABSENT, BASE)

apply_style()
ANA = os.path.join(BASE, "validation", "analysis")
FIG = os.path.join(BASE, "figures", "main_text")

rt = pd.read_csv(os.path.join(BASE, "pangenome/direction_matrices/relay_trait_matrix.tsv"), sep="\t")
LATIN = {"丁酸梭菌": "C. butyricum", "大肠杆菌": "E. coli",
         "加氏乳杆菌": "L. gasseri", "嗜酸乳杆菌": "L. acidophilus",
         "约氏乳杆菌": "L. johnsonii", "双歧杆菌": "Bifidobacterium sp.",
         "罗伊氏乳杆菌": "L. reuteri", "长双歧杆菌": "B. longum",
         "青春双歧杆菌": "B. adolescentis", "普拉梭菌": "F. prausnitzii"}

PATHWAYS = {
    "Lactate-SCFA":    ("LDH/MDH（产乳酸）", "buk（产丁酸）"),
    "Bile acid":       ("BSH（胆盐水解）",   "bai 操纵子（7α-脱羟）"),
    "TMA/TMAO":        ("CutC（产 TMA）",    "TorA（TMAO 呼吸）"),
    "GABA-glutamate":  ("GadA/B（产 GABA）", "GabT（消耗 GABA）"),
}

# relay_trait_matrix holds per-unit CARRIER FRACTIONS (0.12, 0.46, 0.94, ...),
# not 0/1 calls. Positivity is the project-wide threshold of >=60% of a unit's
# genomes; astype(bool) would wrongly score every fractional carrier.
POS = 0.60

fig = plt.figure(figsize=(7.2, 6.2))
gs = fig.add_gridspec(2, 2, hspace=0.56, wspace=0.32)

# ------------------------------------------------ (a) quadrant grid
ax = fig.add_subplot(gs[0, 0])
order = ["U+D+", "U+D−", "U−D+", "U−D−"]
qcol = [NEUTRAL, UP_C, DOWN_C, ABSENT]
grid = np.zeros((len(PATHWAYS), 4))
for i, (pw, (u, d)) in enumerate(PATHWAYS.items()):
    U = rt[u].values >= POS
    D = rt[d].values >= POS
    grid[i] = [(U & D).sum(), (U & ~D).sum(), (~U & D).sum(), (~U & ~D).sum()]
# colour follows the quadrant (categorical), magnitude is printed in the cell
for i in range(len(PATHWAYS)):
    for j in range(4):
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                   facecolor=qcol[j], edgecolor="white", lw=1.2))
        v = int(grid[i, j])
        ax.text(j, i, str(v), ha="center", va="center", fontsize=9,
                fontweight="bold",
                color="white" if j in (1, 2) else "#333333")
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(len(PATHWAYS) - 0.5, -0.5)
ax.set_xticks(range(4))
ax.set_xticklabels([f"{k}\nn = {int(grid[:, j].sum())}" for j, k in enumerate(order)],
                   fontsize=7)
ax.set_yticks(range(len(PATHWAYS)))
ax.set_yticklabels(list(PATHWAYS.keys()), fontsize=7.5)
ax.set_title("(a) Relay configuration by pathway", loc="left")
ax.set_xlabel("unit count")
tot_up_only = int(grid[:, 1].sum())
tot_dn_only = int(grid[:, 2].sum())
ax.text(0.5, -0.30,
        f"upstream-only {tot_up_only} vs downstream-only {tot_dn_only} combinations",
        transform=ax.transAxes, ha="center", va="top", fontsize=7,
        bbox=dict(fc="white", ec=UP_C, lw=0.8, boxstyle="round,pad=0.3"))

# ------------------------------------------------ (b) orphan carriers
ax = fig.add_subplot(gs[0, 1])
orphan_rows = []
for pw, (u, d) in PATHWAYS.items():
    U = rt[u].values >= POS
    D = rt[d].values >= POS
    dn_pos = int(D.sum())
    complete = int((U & D).sum())
    orphan = dn_pos - complete
    orphan_rows.append((pw, dn_pos, complete, orphan,
                        orphan / dn_pos * 100 if dn_pos else np.nan))
xs = np.arange(len(orphan_rows))
carriers = [r[1] for r in orphan_rows]
orph = [r[3] for r in orphan_rows]
ax.bar(xs, carriers, width=0.5, color=NEUTRAL, label="carries downstream")
ax.bar(xs, orph, width=0.5, color=DOWN_C, label="needs external intermediate")
for i, r in enumerate(orphan_rows):
    ax.text(i, r[1] + 0.7, f"{r[3]}/{r[1]}", ha="center", va="bottom",
            fontsize=8, fontweight="bold",
            color=DOWN_C if r[4] >= 50 else "#333333")
    if r[3] > 0:
        ax.text(i, r[3] / 2, f"{r[4]:.0f}%", ha="center", va="center",
                fontsize=8, fontweight="bold", color="white")
ax.set_xticks(xs)
ax.set_xticklabels([r[0] for r in orphan_rows], fontsize=7)
ax.set_ylabel("Units carrying the downstream enzyme")
ax.set_title("(b) Downstream carriers depend on partners", loc="left")
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.set_ylim(0, max(carriers) * 1.42)

# ------------------------------------------------ (c) minimum combination
ax = fig.add_subplot(gs[1, 0])
with open(os.path.join(ANA, "fig7_minimal_combination.json"), encoding="utf-8") as f:
    comb = json.load(f)
single_max = comb["max_single_unit_coverage"]
n_pw = len(PATHWAYS)
# a real comparison: what the best single strain covers vs what the minimal
# combination covers. A single-valued all-green bar would carry no information.
best_names = ["best single\nstrain", f"{len(comb['chosen_units'])}-strain\ncombination"]
best_vals = [single_max, int(comb["pathways_covered"].__len__())]
bars = ax.bar(np.arange(2), best_vals, width=0.46,
              color=[NEUTRAL, BOND_C])
for rect, v in zip(bars, best_vals):
    ax.text(rect.get_x() + rect.get_width() / 2, v + 0.12,
            f"{v}/{n_pw} pathways", ha="center", va="bottom",
            fontsize=8.5, fontweight="bold")
ax.set_xticks(np.arange(2))
ax.set_xticklabels(best_names, fontsize=7.5)
ax.set_ylabel("Pathways with both halves present")
ax.set_ylim(0, n_pw * 1.34)
ax.set_yticks(range(n_pw + 1))
ax.set_title("(c) Multi-strain formulas are required", loc="left")
names = ", ".join(LATIN.get(n, n) for n in comb["chosen_names"])
ax.text(0.5, 0.965,
        f"minimum combination: {names}",
        transform=ax.transAxes, ha="center", va="top", fontsize=7,
        bbox=dict(fc="white", ec=BOND_C, lw=0.8, boxstyle="round,pad=0.3"))
ax.text(0.5, -0.235,
        f"single strains complete at most {single_max} of {n_pw} pathways; "
        f"mean {comb['mean_single_unit_coverage']:.2f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=6.8, color=NEUTRAL)

# ------------------------------------------------ (d) per-unit relay load
ax = fig.add_subplot(gs[1, 1])
n_complete = (((rt[u].values >= POS) & (rt[d].values >= POS))
              for u, d in PATHWAYS.values())
comp_counts = np.sum([x.astype(int) for x in n_complete], axis=0)
vals, cnts = np.unique(comp_counts, return_counts=True)
ax.bar(vals, cnts, width=0.52, color=PALETTE[6])
for v, c in zip(vals, cnts):
    ax.text(v, c + 0.55, str(c), ha="center", va="bottom", fontsize=8,
            fontweight="bold")
ax.set_xlabel("Pathways completed by a single strain")
ax.set_ylabel("Units")
ax.set_xticks(range(len(PATHWAYS)))
ax.set_ylim(0, cnts.max() * 1.34)
ax.set_title("(d) Single strains carry few complete relays", loc="left")
n_zero = int(cnts[vals == 0][0]) if (vals == 0).any() else 0
ax.text(0.5, 0.965,
        f"{n_zero}/{len(rt)} units complete none of the {n_pw} relays",
        transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
        fontweight="bold", color=PALETTE[6],
        bbox=dict(fc="white", ec=PALETTE[6], lw=0.8, boxstyle="round,pad=0.32"))

fig.savefig(os.path.join(FIG, "fig7_first_runner_consequences.png"), dpi=300)
fig.savefig(os.path.join(FIG, "fig7_first_runner_consequences.pdf"))
plt.close(fig)
print("Fig 7 written.")
print(f"  orphan carriers: {[(r[0], f'{r[3]}/{r[1]}') for r in orphan_rows]}")
print(f"  best single strain covers {single_max}/{n_pw}; combination "
      f"{comb['chosen_names']} covers {best_vals[1]}/{n_pw}")
print(f"  {n_zero}/{len(rt)} units complete none")