#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 5 — Clinical corroboration: intermediate accumulation without end-product
         formation  (Plan C, real computation)
=====================================================================
Replaces the old fig11_disrupt_pattern, whose 5x5 disease-by-pathway matrix
was hard-coded (RA / T2D / CD / UC / ASD) and is not derivable from any
dataset in this project. It is replaced with the computation the relay model
actually predicts and the IBDMDB metabolome can test:

    if a relay is broken downstream, then substrates and intermediates should
    accumulate while end products should not be formed.

Every number below is read from relay_metabolite_disease_stats_ibdmdb.tsv.
"""
import sys, os
sys.path.insert(0, r"D:\Probiotics_Pangenome_v5\validation")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from scipy.stats import binomtest
from paper_palette import (apply_style, PALETTE, UP_C, DOWN_C,
                           DISEASE_COLORS, NEUTRAL, ABSENT, BASE)

# clinical group colours are fixed by the global palette so that CD and UC keep
# the same colour in every figure of the paper
DIS_C = DISEASE_COLORS["CD"]     # slot 2  red
CTRL_C = DISEASE_COLORS["UC"]    # slot 12 salmon

apply_style()
ANA = os.path.join(BASE, "validation", "analysis")
FIG = os.path.join(BASE, "figures", "main_text")

m = pd.read_csv(os.path.join(ANA, "relay_metabolite_disease_stats_ibdmdb.tsv"), sep="\t")
m["lfc_cd"] = m.apply(lambda r: np.log2((r.CD_mean + 1) / (r.nonIBD_mean + 1)), axis=1)
m["lfc_uc"] = m.apply(lambda r: np.log2((r.UC_mean + 1) / (r.nonIBD_mean + 1)), axis=1)
m["pmin"] = m[["p_CD_vs_nonIBD", "p_UC_vs_nonIBD"]].min(axis=1)
sig = m[m.pmin < 0.05].copy().sort_values("lfc_cd", ascending=False)

ROLE_ORDER = ["substrate", "intermediate", "end_product"]
ROLE_LABEL = {"substrate": "Substrate", "intermediate": "Intermediate",
              "end_product": "End product"}
ROLE_COL = {"substrate": PALETTE[8], "intermediate": UP_C, "end_product": DOWN_C}

PW_LABEL = {"bile_acid": "Bile acid", "tma_tmao": "TMA/TMAO",
            "tryptophan": "Tryptophan", "lactate_scfa": "Lactate-SCFA",
            "gaba_glutamate": "GABA-glutamate"}

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.32, 1], hspace=0.72, wspace=0.34)

# ------------------------------------------------ (a) stratified log2FC ladder
ax = fig.add_subplot(gs[0, :])
y = np.arange(len(sig))[::-1]
for yy, (_, r) in zip(y, sig.iterrows()):
    ax.barh(yy, r.lfc_cd, height=0.62,
            color=ROLE_COL[r.role] if r.pmin < 0.05 else ABSENT,
            edgecolor="none")
    # UC value as a tick, so both diseases are visible per metabolite
    ax.plot([r.lfc_uc], [yy], marker="|", ms=7, mew=1.3,
            color=CTRL_C, zorder=5)
ax.axvline(0, color="#555555", lw=0.9)
ax.set_yticks(y)
ax.set_yticklabels([f"{PW_LABEL[r.pathway]}: {r.metabolite}" for _, r in sig.iterrows()],
                   fontsize=6.5)
ax.set_xlabel("log$_2$ fold change in disease vs non-IBD")
ax.set_title("(a) Relay-step-resolved metabolite change (significant metabolites, "
             f"n = {len(sig)} of {len(m)})", loc="left")
ax.set_xlim(-2.4, 3.9)
ax.tick_params(axis="x", labelsize=7)
# group separators + role brackets
prev_role, ystart = None, None
for yy, (_, r) in zip(y, sig.iterrows()):
    if r.role != prev_role:
        if prev_role is not None:
            ax.axhline(ystart - 0.5, color=NEUTRAL, lw=0.6, ls=":")
        prev_role, ystart = r.role, yy
ax.axhline(ystart - 0.5, color=NEUTRAL, lw=0.6, ls=":")
handles = [Patch(fc=ROLE_COL[r], ec="none", label=ROLE_LABEL[r]) for r in ROLE_ORDER]
handles.append(plt.Line2D([], [], marker="|", ls="none", color=CTRL_C,
                          mew=1.6, ms=9, label="Ulcerative colitis"))
handles.append(Patch(fc=ABSENT, ec="none", label="n.s."))
ax.legend(handles=handles, frameon=False, ncol=5, loc="lower right",
          bbox_to_anchor=(1.0, -0.30), handlelength=1.0, columnspacing=1.3)

# ------------------------------------------------ (b) stratified medians
ax = fig.add_subplot(gs[1, 0])
stats_rows = []
for role in ROLE_ORDER:
    s_all = m[m.role == role]
    s_sig = sig[sig.role == role]
    stats_rows.append((role, len(s_sig), len(s_all),
                       s_sig.lfc_cd.median(), s_sig.lfc_uc.median()))
medians = [r[3] for r in stats_rows]
x = np.arange(len(ROLE_ORDER))
w = 0.36
b1 = ax.bar(x - w / 2, medians, w, color=DIS_C, label="Crohn's disease")
b2 = ax.bar(x + w / 2, [r[4] for r in stats_rows], w, color=CTRL_C,
            label="Ulcerative colitis")
ax.axhline(0, color="#555555", lw=0.9)
ax.set_xticks(x)
ax.set_xticklabels([f"{ROLE_LABEL[r]}\n({s}/{n} sig.)" for r, s, n, _, _ in stats_rows],
                   fontsize=7)
ax.set_ylabel("Median log$_2$FC vs non-IBD")
ax.set_title("(b) Change accumulates upstream of the block", loc="left")
ax.legend(frameon=False, fontsize=7, loc="upper right")
ax.set_ylim(min(-0.6, min(medians) * 1.5), max(medians) * 1.55 + 0.5)
for rect, v in list(zip(b1, medians)) + list(zip(b2, [r[4] for r in stats_rows])):
    off = 0.08 if v >= 0 else -0.08
    ax.text(rect.get_x() + rect.get_width() / 2, v + off, f"{v:+.2f}",
            ha="center", va="bottom" if v >= 0 else "top", fontsize=6.5)

si = sig[sig.role.isin(["substrate", "intermediate"])]
ep = sig[sig.role == "end_product"]
mw_p = stats.mannwhitneyu(si.lfc_cd.abs(), ep.lfc_cd.abs()).pvalue
cd_conc = int(((sig.role.isin(["substrate", "intermediate"])) & (sig.lfc_cd > 0)).sum()
              + ((sig.role == "end_product") & (sig.lfc_cd < 0)).sum())
cd_bt = binomtest(cd_conc, len(sig), 0.5, alternative="greater").pvalue
ax.text(0.5, -0.235,
        f"|log$_2$FC|: substrate+intermediate {si.lfc_cd.abs().median():.2f} "
        f"vs end product {ep.lfc_cd.abs().median():.2f}   MWU p = {mw_p:.3f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=6.8,
        bbox=dict(fc="white", ec=NEUTRAL, lw=0.6, boxstyle="round,pad=0.3"))

# ------------------------------------------------ (c) directional concordance
ax = fig.add_subplot(gs[1, 1])
cats = ["CD", "UC"]
vals = []
for col in ["lfc_cd", "lfc_uc"]:
    v = int(((sig.role.isin(["substrate", "intermediate"])) & (sig[col] > 0)).sum()
            + ((sig.role == "end_product") & (sig[col] < 0)).sum())
    vals.append(v)
n = len(sig)
xs = np.arange(2)
bars = ax.bar(xs, vals, width=0.44, color=[DIS_C, CTRL_C])
for rect, v in zip(bars, vals):
    ax.text(rect.get_x() + rect.get_width() / 2, v + 0.7,
            f"{v}/{n}\n({v/n*100:.0f}%)", ha="center", va="bottom", fontsize=8,
            fontweight="bold")
ax.axhline(n / 2, color=NEUTRAL, lw=1.0, ls="--")
ax.text(-0.44, n / 2, "chance", fontsize=6.5, color=NEUTRAL,
        va="bottom", ha="left")
ax.set_xticks(xs)
ax.set_xticklabels(["Crohn's disease", "Ulcerative colitis"], fontsize=7.5)
ax.set_ylabel("Metabolites with relay-predicted direction")
ax.set_ylim(0, n * 1.34)
ax.set_yticks(np.arange(0, n + 1, 4))
ax.set_title("(c) Direction matches the relay prediction", loc="left")
bt_p = binomtest(vals[0], n, 0.5, alternative="greater").pvalue
ax.text(0.5, -0.235,
        f"binomial test vs chance: CD p = {cd_bt:.4f}, UC p = "
        f"{binomtest(vals[1], n, 0.5, alternative='greater').pvalue:.4f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=6.8,
        bbox=dict(fc="white", ec=NEUTRAL, lw=0.6, boxstyle="round,pad=0.3"))

fig.savefig(os.path.join(FIG, "fig5_relay_disruption_clinical.png"), dpi=300)
fig.savefig(os.path.join(FIG, "fig5_relay_disruption_clinical.pdf"))
plt.close(fig)

with open(os.path.join(ANA, "fig5_relay_disruption_stats.txt"), "w", encoding="utf-8") as f:
    f.write("Fig5 relay disruption stats (IBDMDB, 21 metagenomes/metabolomes)\n")
    f.write(f"significant metabolites (p<0.05 vs nonIBD in either disease): "
            f"{len(sig)}/{len(m)}\n")
    for role, s, nn, mc, mu in stats_rows:
        f.write(f"  {role:13s} sig {s}/{nn}  medianCD={mc:+.3f}  medianUC={mu:+.3f}\n")
    f.write(f"|log2FC| substrate+intermediate {si.lfc_cd.abs().median():.3f} "
            f"(n={len(si)}) vs end_product {ep.lfc_cd.abs().median():.3f} "
            f"(n={len(ep)}); MWU p={mw_p:.4g}\n")
    f.write(f"concordant direction: CD {vals[0]}/{n} p={cd_bt:.4g}; "
            f"UC {vals[1]}/{n} p={binomtest(vals[1], n, 0.5, alternative='greater').pvalue:.4g}\n")
    ba = sig[sig.pathway == "bile_acid"]
    f.write(f"bile acid: primary/intermediate median {ba[ba.role=='intermediate'].lfc_cd.median():+.3f}"
            f" (n={(ba.role=='intermediate').sum()}) vs secondary/end "
            f"{ba[ba.role=='end_product'].lfc_cd.median():+.3f} "
            f"(n={(ba.role=='end_product').sum()}); "
            f"MWU p={stats.mannwhitneyu(ba[ba.role=='intermediate'].lfc_cd, ba[ba.role=='end_product'].lfc_cd).pvalue:.4g}\n")
print("Fig 5 written.")
print(f"  sig {len(sig)}/{len(m)}  MWU p={mw_p:.4f}  CD {vals[0]}/{n} p={cd_bt:.4g}"
      f"  UC {vals[1]}/{n}")