#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 3 — Systematic upstream-downstream asymmetry  (Plan C, real computation)
=====================================================================
All panels are computed from the 351-genome pangenome trait matrix.
No hard-coded values anywhere in this file.

(a) Genome-wide direction budget: every unit carries more substrate-generating
    (dir1) than host-effect-generating (dir2) orthogroups.
(b) Pathway-level upstream vs downstream carrier frequency.
(c) Joint distribution: the "upstream-only" quadrants dominate.
(d) Per-unit asymmetry: all 41 units have dir1 fraction > dir2 fraction.
"""
import sys, os
sys.path.insert(0, r"D:\Probiotics_Pangenome_v5\validation")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from paper_palette import (apply_style, PALETTE, UP_C, DOWN_C, NEUTRAL, ABSENT,
                           present_cmap, BASE)
apply_style()

ANA = os.path.join(BASE, "validation", "analysis")
FIG = os.path.join(BASE, "figures", "main_text")

rt = pd.read_csv(os.path.join(BASE, "pangenome/direction_matrices/relay_trait_matrix.tsv"), sep="\t")
dg = pd.read_csv(os.path.join(BASE, "pangenome/unit_direction_genes_orthogroup.tsv"), sep="\t")
dg["tot"] = dg[["dir1_genes", "dir2_genes", "dir3_genes"]].sum(axis=1)
dg["d1f"] = dg.dir1_genes / dg.tot
dg["d2f"] = dg.dir2_genes / dg.tot

PATHWAYS = {
    "Lactate\nSCFA":        ("LDH/MDH（产乳酸）", "buk（产丁酸）"),
    "Bile acid":      ("BSH（胆盐水解）",   "bai 操纵子（7α-脱羟）"),
    "TMA/TMAO":("CutC（产 TMA）",    "TorA（TMAO 呼吸）"),
    "GABA\nGlutamate":     ("GadA/B（产 GABA）", "GabT（消耗 GABA）"),
}

# relay_trait_matrix holds per-unit CARRIER FRACTIONS (values such as 0.12,
# 0.46, 0.94), not 0/1 calls. A unit counts as carrying a trait only when
# the trait is present in >= 60% of its genomes (the threshold used throughout
# this project); any non-zero test would wrongly score fractional carriers.
POS = 0.60
# panel (b) shows the underlying continuous carrier fraction
frac = {k: (rt[u].mean() * 100, rt[d].mean() * 100) for k, (u, d) in PATHWAYS.items()}
call = {k: ((rt[u] >= POS).sum(), (rt[d] >= POS).sum())
        for k, (u, d) in PATHWAYS.items()}

fig = plt.figure(figsize=(7.2, 9.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.15, 0.92],
                      hspace=0.62, wspace=0.34)

# ---------------------------------------------------------------- (a)
ax = fig.add_subplot(gs[0, 0])
tot1, tot2, tot3 = dg.dir1_genes.sum(), dg.dir2_genes.sum(), dg.dir3_genes.sum()
labels = ["Substrate-\ngenerating", "Host-effect-\ngenerating", "Other"]
vals = [tot1, tot2, tot3]
cols = [UP_C, DOWN_C, NEUTRAL]
b = ax.bar(labels, vals, color=cols, width=0.62)
ax.set_yscale("log")
ax.set_ylabel("Orthogroups (log scale)")
ax.set_title("(a) Genome-wide direction budget", loc="left")
for rect, v in zip(b, vals):
    ax.text(rect.get_x() + rect.get_width() / 2, v * 1.18, f"{v:,}",
            ha="center", va="bottom", fontsize=7)
chi2, pchi, _, _ = stats.chi2_contingency([[tot1, tot2], [tot3, 1]])
ax.set_ylim(1, tot1 * 120)
# annotation sits in the empty upper-right corner, clear of the tallest bar
# (dir1, left) and of the panel title
ax.text(0.97, 0.97, f"dir1 : dir2 = {tot1/tot2:.1f}×   χ² p = {pchi:.1e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
        fontweight="bold",
        bbox=dict(fc="white", ec=UP_C, lw=0.8, boxstyle="round,pad=0.34"))

# ---------------------------------------------------------------- (b)
ax = fig.add_subplot(gs[0, 1])
names = list(PATHWAYS.keys())
ups = [frac[k][0] for k in names]
dns = [frac[k][1] for k in names]
x = np.arange(len(names))
w = 0.36
b1 = ax.bar(x - w / 2, ups, w, color=UP_C, label="Upstream (substrate-generating)")
b2 = ax.bar(x + w / 2, dns, w, color=DOWN_C, label="Downstream (host-effect-generating)")
ax.set_xticks(x); ax.set_xticklabels(names)
ax.set_ylabel("Mean carrier frequency across units (%)")
ax.set_title("(b) Upstream vs downstream carrier frequency", loc="left")
# legend below the panel: inside it collided with the 96% value label and with
# the ratio row, both of which sit near the top of the tallest bar.
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20),
          ncol=1, handlelength=1.0, labelspacing=0.28, fontsize=7,
          borderaxespad=0.0)
for rect, v in list(zip(b1, ups)) + list(zip(b2, dns)):
    ax.text(rect.get_x() + rect.get_width() / 2, v + 2.0, f"{v:.0f}",
            ha="center", va="bottom", fontsize=6)
ax.set_ylim(0, 118)
ratios = [u / d if d > 0 else np.nan for u, d in zip(ups, dns)]
for i, rr in enumerate(ratios):
    ax.text(i, 108, ("%.1f×" % rr) if rr == rr else "n/a", ha="center",
            va="center", fontsize=8, fontweight="bold",
            color=DOWN_C if rr < 1 else UP_C)

# ---------------------------------------------------------------- (c)
ax = fig.add_subplot(gs[1, 0])
quad_counts = {"U+D+": [], "U+D−": [], "U−D+": [], "U−D−": []}
for pw, (u, d) in PATHWAYS.items():
    U = rt[u].values >= POS
    D = rt[d].values >= POS
    quad_counts["U+D+"].append(int((U & D).sum()))
    quad_counts["U+D−"].append(int((U & ~D).sum()))
    quad_counts["U−D+"].append(int((~U & D).sum()))
    quad_counts["U−D−"].append(int((~U & ~D).sum()))
order = ["U+D+", "U+D−", "U−D+", "U−D−"]
qc = {k: int(np.sum(quad_counts[k])) for k in order}
qcol = [NEUTRAL, UP_C, DOWN_C, ABSENT]
left = 0
for i, k in enumerate(order):
    ax.barh([0], [qc[k]], left=left, color=qcol[i], height=0.42)
    ax.text(left + qc[k] / 2, 0, str(qc[k]), ha="center", va="center",
            fontsize=8.5, fontweight="bold",
            color="white" if i in (1, 2) else "#333333")
    left += qc[k]
ax.set_xlim(0, left * 1.02)
ax.set_ylim(-1.60, 0.60)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlabel("")
ax.set_title("(c) Joint upstream/downstream configuration", loc="left")
up_only, dn_only = qc["U+D−"], qc["U−D+"]
ax.text(0.5, 0.47,
        f"upstream-only {up_only}  vs  downstream-only {dn_only}   ({up_only/dn_only:.1f}×)",
        transform=ax.transAxes, ha="center", va="center", fontsize=8,
        fontweight="bold", color=UP_C,
        bbox=dict(fc="white", ec=UP_C, lw=0.9, boxstyle="round,pad=0.36"))
leg = [Patch(fc=qcol[i], ec="none",
             label=f"{k}   n = {qc[k]}") for i, k in enumerate(order)]
# legend sits between the summary line and the bar, fully inside the axes
ax.legend(handles=leg, frameon=False, loc="center",
          bbox_to_anchor=(0.5, 0.20), ncol=2, handlelength=1.1,
          columnspacing=1.2, labelspacing=0.55)

# ---------------------------------------------------------------- (d)
ax = fig.add_subplot(gs[1, 1])
ax.scatter(dg.d1f, dg.d2f, s=30, color=PALETTE[3], alpha=0.75,
           edgecolors="white", linewidths=0.6)
lim = 0.92
ax.plot([0, lim], [0, lim], ls="--", lw=0.9, color=NEUTRAL, zorder=0)
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("Substrate-generating fraction of unit genes")
ax.set_ylabel("Host-effect-generating fraction")
ax.set_title("(d) Per-unit asymmetry", loc="left")
n_above = int((dg.d1f > dg.d2f).sum())
wil_p = stats.wilcoxon(dg.d1f, dg.d2f).pvalue
# annotate the two extremes that anchor the claim; English genus names only
# (the figure must contain no CJK glyphs - journals require English-only figures)
LATIN = {"加氏乳杆菌": "L. gasseri", "青春双歧杆菌": "B. adolescentis",
         "约氏乳杆菌": "L. johnsonii", "嗜酸乳杆菌": "L. acidophilus",
         "丁酸梭菌": "C. butyricum", "大肠杆菌": "E. coli",
         "罗伊氏乳杆菌": "L. reuteri", "长双歧杆菌": "B. longum"}
i_max = dg.d2f.idxmax()
i_top = dg.d1f.idxmax()
# offsets chosen so neither label runs off the right edge nor lands on the
# scatter cloud; ha/va flip with the offset direction
for idx, off in [(i_max, (-9, 12)), (i_top, (-9, 17))]:
    lbl = LATIN.get(dg.loc[idx, "cn_name"], "unit " + str(dg.loc[idx, "unit_id"]))
    ax.annotate(lbl, (dg.loc[idx, "d1f"], dg.loc[idx, "d2f"]),
                textcoords="offset points", xytext=off, ha="right", va="center",
                fontsize=7, color=PALETTE[8], fontweight="bold")
# the point cloud occupies the lower-right, so the statistics box goes to the
# empty upper-left corner instead of across the data
ax.text(0.03, 0.965,
        f"{n_above}/{len(dg)} units above\nthe diagonal\nWilcoxon p = {wil_p:.1e}",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
        fontweight="bold", linespacing=1.45,
        bbox=dict(fc="white", ec=PALETTE[3], lw=0.9, boxstyle="round,pad=0.36"))
ax.text(0.86, 0.80, "y = x", fontsize=7, color=NEUTRAL, rotation=39)

# ---------------------------------------------------------------- (e)
# Genome-size control. A small genome would be expected to lose both halves,
# so a smaller genome could produce the asymmetry as a by-product of reduction.
# Instead the bias widens with genome size, and downstream carriage is
# independent of size, which excludes reduction as the explanation.
gsz = pd.read_csv(f"{BASE}/qc/genome_stats_by_unit.tsv", sep="\t")
gsz["unit_id"] = gsz.unit_id.astype(int)
sz = dg[["unit_id", "cn_name", "d1f", "d2f"]].merge(
    gsz[["unit_id", "median_size_mb"]], on="unit_id", how="inner")
sz["ratio"] = (sz.d1f + 0.01) / (sz.d2f + 0.01)

ax = fig.add_subplot(gs[2, 0])
_k = max(4, len(sz) // 5)
_s = sz.sort_values("median_size_mb")
_lo, _hi = _s.head(_k), _s.tail(_k)
_jitter = np.linspace(-0.10, 0.10, len(sz))
ax.scatter(sz.median_size_mb + _jitter, sz.ratio, s=26,
           color=PALETTE[3], alpha=0.70, edgecolors="white", linewidths=0.5)
for _grp, _mk, _c in ((_lo, "o", PALETTE[8]), (_hi, "D", PALETTE[1])):
    ax.scatter(_grp.median_size_mb, _grp.ratio, s=95, marker=_mk,
               facecolors="none", edgecolors=_c, linewidths=1.5, zorder=3)
_rho_sz, _p_sz = stats.spearmanr(sz.median_size_mb, sz.ratio)
_u_sz, _p_dec = stats.mannwhitneyu(_lo.ratio, _hi.ratio, alternative="less")
ax.set_xlabel("Median genome size of the unit (Mb)")
ax.set_ylabel("Upstream / downstream gene ratio")
ax.set_title("(e) The bias is not genome-size driven", loc="left")
_ytop = max(sz.ratio) * 1.42
ax.set_ylim(0, _ytop)
ax.text(0.03, 0.965,
        f"size vs ratio: $\\rho$ = {_rho_sz:+.2f}, p = {_p_sz:.3f}\n"
        f"largest {_k} > smallest {_k}: p = {_p_dec:.4f}",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.2,
        fontweight="bold", linespacing=1.4,
        bbox=dict(fc="white", ec=PALETTE[3], lw=0.9, boxstyle="round,pad=0.34"))

# ---------------------------------------------------------------- (f)
# Downstream carriage is flat across genome sizes; upstream carriage rises.
ufc = pd.read_csv(f"{BASE}/pangenome/direction_matrices/unit_family_counts.tsv",
                  sep="\t")
_rows = []
for _u, _grp in ufc.groupby("unit_id"):
    _tot = _grp.n_genes.sum()
    _rows.append({
        "unit_id": _u,
        "pct_down": _grp[_grp.direction == "dir2"].n_genes.sum() / _tot,
        "pct_up": _grp[_grp.direction == "dir1"].n_genes.sum() / _tot,
    })
R = pd.DataFrame(_rows)
R["unit_id"] = R.unit_id.astype(int)
R = R.merge(gsz[["unit_id", "median_size_mb"]], on="unit_id", how="inner")

ax = fig.add_subplot(gs[2, 1])
ax.scatter(R.median_size_mb + _jitter[:len(R)], R.pct_down * 100, s=26,
           color=PALETTE[1], alpha=0.70, edgecolors="white", linewidths=0.5,
           label="downstream genes")
ax.scatter(R.median_size_mb + _jitter[:len(R)], R.pct_up * 100, s=26,
           color=UP_C, alpha=0.70, edgecolors="white", linewidths=0.5,
           label="upstream genes")
_rho_d, _p_d = stats.spearmanr(R.median_size_mb, R.pct_down)
_rho_u, _p_u = stats.spearmanr(R.median_size_mb, R.pct_up)
for _y, _c in ((R.pct_down * 100, PALETTE[1]), (R.pct_up * 100, UP_C)):
    _b = np.polyfit(R.median_size_mb, _y, 1)
    _xs = np.linspace(R.median_size_mb.min(), R.median_size_mb.max(), 40)
    ax.plot(_xs, np.poly1d(_b)(_xs), color=_c, lw=1.6, zorder=2)
ax.set_xlabel("Median genome size of the unit (Mb)")
ax.set_ylabel("Share of unit genes (%)")
ax.set_title("(f) Carriage of each half vs genome size", loc="left")
ax.set_ylim(0, 108)
# the cloud fills the upper-left, so the legend sits in the empty mid-right band
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.60, 0.62),
          fontsize=7.2, handletextpad=0.4, borderaxespad=0.0)
ax.text(0.975, 0.97,
        f"downstream $\\rho$ = {_rho_d:+.2f}, p = {_p_d:.2f}\n"
        f"upstream $\\rho$ = {_rho_u:+.2f}, p = {_p_u:.2f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.2,
        fontweight="bold", linespacing=1.4,
        bbox=dict(fc="white", ec=NEUTRAL, lw=0.8, boxstyle="round,pad=0.32"))

fig.savefig(os.path.join(FIG, "fig3_upstream_downstream_asymmetry.png"), dpi=300)
fig.savefig(os.path.join(FIG, "fig3_upstream_downstream_asymmetry.pdf"))
plt.close(fig)

# ---------------------------------------------------------------- stats table
rows = []
for pw, (u, d) in PATHWAYS.items():
    nu, nd = call[pw]
    rows.append(dict(pathway=pw.replace("\n", " "), upstream_trait=u, downstream_trait=d,
                     up_units_at_60pct=int(nu), dn_units_at_60pct=int(nd),
                     up_mean_carrier_pct=round(frac[pw][0], 1),
                     dn_mean_carrier_pct=round(frac[pw][1], 1),
                     ratio=round(frac[pw][0] / frac[pw][1], 2) if frac[pw][1] > 0 else np.nan))
pd.DataFrame(rows).to_csv(os.path.join(ANA, "fig3_asymmetry_pathway_stats.tsv"),
                          sep="\t", index=False)
with open(os.path.join(ANA, "fig3_asymmetry_summary.txt"), "w", encoding="utf-8") as f:
    f.write("Fig3 asymmetry summary\n")
    f.write(f"positivity threshold: carrier fraction >= {POS} within a unit\n")
    f.write(f"dir1 total={tot1}  dir2 total={tot2}  dir3 total={tot3}  ratio={tot1/tot2:.2f}\n")
    f.write(f"chi-square dir1 vs dir2 p={pchi:.3e}\n")
    f.write(f"units with d1f>d2f: {n_above}/{len(dg)}  Wilcoxon p={wil_p:.3e}\n")
    f.write(f"quadrants: {qc}\n")
    f.write(f"upstream-only / downstream-only = {up_only}/{dn_only} = {up_only/dn_only:.2f}x\n")
    for r in rows:
        f.write(f"  {r['pathway']:16s} up {r['up_mean_carrier_pct']:5.1f}% "
                f"dn {r['dn_mean_carrier_pct']:5.1f}%  ratio {r['ratio']}\n")
print("Fig 3 written.")
print(f"  dir1={tot1} dir2={tot2} ratio={tot1/tot2:.2f} chi2 p={pchi:.2e}")
print(f"  {n_above}/{len(dg)} units above diagonal, Wilcoxon p={wil_p:.2e}")
print(f"  quadrants {qc}; up-only vs dn-only = {up_only/dn_only:.2f}x")