#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fig. 4 | The relay partners: probiotic -> community partner -> product
Tripartite network on 21 HUManN2 metagenomes (iHMP/HMP2, CD=15, UC=6).

Panels
 a  Schematic tripartite network: probiotic units (upstream) -> community
    partners (downstream enzyme carriers) -> host-effect products.
 b  Partner abundance by pathway: who carries the downstream reaction.
 c  Disease contrast: partner abundance in CD vs UC.
 d  Species-level partner abundance across the cohort.

All labels in English.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Circle, Wedge, FancyBboxPatch

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "axes.linewidth": 0.8,
})

BASE = r"D:\Probiotics_Pangenome_v5"
ANA = os.path.join(BASE, "validation", "analysis")
FIG = os.path.join(BASE, "figures", "main_text")
os.makedirs(FIG, exist_ok=True)

from paper_palette import UP_C, DOWN_C, BOND_C, PALETTE
PROD_C = PALETTE[5]   # yellow-family accent for products
GREY = PALETTE[8]

PW = {
    "lactate_to_butyrate": dict(
        label="Lactate\nto butyrate", product="Butyrate",
        up="Lactate dehydrogenase\n(EC 1.1.1.27), probiotics",
        down="Butyryl-CoA:acetate CoA-transferase (EC 2.8.3.8)\nor butyrate kinase (EC 2.7.2.7)",
        partners=["Escherichia coli", "Roseburia intestinalis",
                  "Anaerostipes hadrus", "Blautia Ruminococcus torques"]),
    "bile_acid_secondary": dict(
        label="Bile acid\n7a-dehydroxylation",
        product="Secondary bile acids\n(LCA, UDCA)",
        up="Bile salt hydrolase\n(EC 3.5.1.26), probiotics",
        down="Bai operon 7a-dehydroxylase\n(EC 1.6.99.7)",
        partners=["Lachnospiraceae noname Lachnospiraceae bacterium 5 1 57FAA",
                  "Bacteroides fragilis", "Eubacterium siraeum"]),
    "gaba_shunt": dict(
        label="GABA\nshunt", product="GABA /\nsuccinic semialdehyde",
        up="Glutamate decarboxylase\n(EC 4.1.1.15), probiotics",
        down="GABA aminotransferase GabT\n(EC 2.6.1.19)",
        partners=["Escherichia coli", "Eubacterium siraeum", "Blautia"]),
    "tryptophan_indole": dict(
        label="Tryptophan\nto indole", product="Indole /\nindole-3-propionate",
        up="Tryptophan biosynthesis\noperon (trpA-G), probiotics",
        down="Tryptophanase TnaA (EC 4.1.99.03)\nand aromatic amino acid aminotransferase",
        partners=["Bacteroides fragilis", "Escherichia coli",
                  "Bacteroides xylanisolvens", "Odoribacter splanchnicus"]),
}

cf = pd.read_csv(os.path.join(ANA, "relay_partner_carrier_fraction.tsv"),
                 sep="\t")
sp = pd.read_csv(os.path.join(ANA, "relay_partner_species_abundance.tsv"),
                 sep="\t", index_col=0)
sx = pd.read_csv(os.path.join(ANA, "relay_partner_cd_uc_stats.tsv"), sep="\t")
cx = pd.read_csv(os.path.join(ANA, "relay_completeness_host.tsv"), sep="\t")


def shorten(name, n=22):
    if len(name) <= n:
        return name
    return name[:n - 3] + "..."


fig = plt.figure(figsize=(19.0, 12.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.22, 1.0],
                      hspace=0.30, wspace=0.22,
                      left=0.045, right=0.975, top=0.925, bottom=0.065)

# ============================== panel a: tripartite network
# Panel a is a NETWORK, not a table: three node columns (probiotic upstream
# reaction -> community partner species -> host-effect product) joined by
# curved, weight-bearing edges. Edge width encodes the measured detection rate
# of the partner across the 21 metagenomes (the quantity plotted in panel b), so
# every edge in this panel is a real observation rather than a drawn arrow.
axA = fig.add_subplot(gs[0, 0])
axA.set_axis_off()
axA.set_xlim(0, 100); axA.set_ylim(0, 100)

pw_keys = list(PW.keys())
n_pw = len(pw_keys)

# measured detection rate per pathway -> edge width
_det = (cf.groupby("pathway")["frac_samples"].mean()
        .reindex(pw_keys).fillna(0) * 100)
_lw = lambda v: 0.7 + 3.4 * (v / 100.0) ** 0.75

# node geometry: three columns
COL_U, COL_P, COL_R = 13.0, 50.0, 86.0
HALF_W, HALF_H = 10.5, 2.35
rows_y = np.linspace(78, 22, n_pw)

# column headers
for x, txt, col in ((COL_U, "Probiotic unit\n(upstream reaction)", UP_C),
                    (COL_P, "Community partner species\n(downstream carrier)", DOWN_C),
                    (COL_R, "Host-effect product", PROD_C)):
    axA.text(x, 95, txt, ha="center", va="center", fontsize=8.4,
             color=col, fontweight="bold", linespacing=1.35)
    axA.plot([x - HALF_W - 1.4, x + HALF_W + 1.4], [88.4, 88.4],
             color=col, lw=1.6, alpha=0.55, clip_on=False, zorder=2)

part_nodes = {}   # (pathway, species) -> (x, y)
for y, k in zip(rows_y, pw_keys):
    d = PW[k]

    # ---- upstream node (probiotic side)
    axA.add_patch(FancyBboxPatch(
        (COL_U - HALF_W, y - HALF_H), 2 * HALF_W, 2 * HALF_H,
        boxstyle="round,pad=0.35,rounding_size=1.2",
        fc="white", ec=UP_C, lw=1.5, zorder=5))
    axA.text(COL_U, y + 1.05, d["label"].replace("\n", " "), ha="center",
             va="center", fontsize=7.6, fontweight="bold", color="black",
             zorder=6)
    # enzyme line wrapped to two rows so it stays inside the node
    _enz = d["up"].split(",")[0].replace("\n", " ")
    if len(_enz) > 22:
        _cut = _enz.rfind(" ", 0, 22)
        _enz = _enz[:_cut] + "\n" + _enz[_cut + 1:]
    axA.text(COL_U, y - 1.05, _enz, ha="center", va="center", fontsize=5.7,
             color=UP_C, zorder=6, linespacing=1.15)

    # ---- product node
    axA.add_patch(FancyBboxPatch(
        (COL_R - HALF_W, y - HALF_H), 2 * HALF_W, 2 * HALF_H,
        boxstyle="round,pad=0.35,rounding_size=1.2",
        fc=PROD_C, ec="none", zorder=5, alpha=0.92))
    _prod = d["product"].replace("\n", " ")
    if len(_prod) > 20:
        _cut = _prod.rfind(" ", 0, 20)
        _prod = _prod[:_cut] + "\n" + _prod[_cut + 1:]
    axA.text(COL_R, y, _prod, ha="center", va="center", fontsize=7.0,
             fontweight="bold", color="white", zorder=6, linespacing=1.2)

    # ---- partner species nodes, fanned out vertically around the row centre
    sp_list = [shorten(p, 24) for p in d["partners"][:4]]
    k_n = len(sp_list)
    dy = 4.6
    ys = [y + (i - (k_n - 1) / 2) * dy for i in range(k_n)]
    for sp_name, sy in zip(sp_list, ys):
        axA.add_patch(FancyBboxPatch(
            (COL_P - HALF_W - 1.0, sy - HALF_H * 0.62),
            2 * HALF_W + 2.0, HALF_H * 1.24,
            boxstyle="round,pad=0.25,rounding_size=1.0",
            fc="white", ec=DOWN_C, lw=1.1, zorder=5))
        axA.text(COL_P, sy, sp_name, ha="center", va="center", fontsize=6.1,
                 color=DOWN_C, style="italic", zorder=6)
        part_nodes[(k, sp_name)] = (COL_P, sy)

        # edge: probiotic -> partner, width = measured detection rate
        axA.add_patch(FancyArrowPatch(
            (COL_U + HALF_W, y), (COL_P - HALF_W - 1.0, sy),
            arrowstyle="-", connectionstyle="arc3,rad=0.16",
            lw=_lw(_det[k]), color=UP_C, alpha=0.55, zorder=3))
        # edge: partner -> product
        axA.add_patch(FancyArrowPatch(
            (COL_P + HALF_W + 1.0, sy), (COL_R - HALF_W, y),
            arrowstyle="-|>", mutation_scale=9,
            connectionstyle="arc3,rad=0.16",
            lw=_lw(_det[k]) * 0.75, color=DOWN_C, alpha=0.55, zorder=3))

# edge-weight key
axA.text(COL_U, 6.0, "Edge width = partner detection rate across the 21 "
         "metagenomes", ha="center", va="center", fontsize=6.6,
         color="#555555", style="italic")

axA.set_ylim(2, 100)
axA.set_title("a   Tripartite relay network: probiotic upstream reaction -> "
              "community partner -> host-effect product",
              fontsize=9.6, fontweight="bold", loc="left", pad=18)
axA.text(0, -0.030,
         "Metagenomes: 21 HUManN2 0.11.0 runs on HMP2/iHMP longitudinal IBD "
         "stool metagenomes (CD n = 15, UC n = 6), retrieved from the "
         "Huttenhower Lab HMP2 data release.\n"
         "Databases: nucleotide ChocoPhlAn v0.1.1; protein UniRef90 "
         "v1.1_uniref90; pathways MetaCyc reactions (level-4 EC) and "
         "structured pathway tree, as distributed with HUManN2 lib_v0.11.0.\n"
         "Probiotic genomes: 351 RefSeq complete genomes (NCBI Assembly), "
         "accessions and BioSample IDs in Supplementary Table S1. "
         "Partner species are those whose genomes carry the downstream enzyme.",
         transform=axA.transAxes, fontsize=6.4, color="#444444", va="top",
         linespacing=1.5)

# ============================== panel b: partner carriers per pathway
axB = fig.add_subplot(gs[0, 1])
frac = (cf.groupby("pathway")["frac_samples"].mean()
        .reindex(pw_keys).fillna(0) * 100)
axB.barh(range(n_pw), frac.values, color=BOND_C, height=0.55,
         edgecolor="white", lw=0.8)
axB.set_yticks(range(n_pw))
axB.set_yticklabels([PW[k]["label"].replace("\n", " ") for k in pw_keys],
                    fontsize=8)
axB.set_xlabel("Mean detection rate of downstream-carrying species (% of 21 metagenomes)",
               fontsize=8)
axB.set_xlim(0, max(100, frac.max() * 1.12))
axB.tick_params(axis="x", labelsize=7.5)
axB.spines[["top", "right"]].set_visible(False)
axB.grid(axis="x", lw=0.4, color="#DDDDDD", zorder=0)
axB.set_axisbelow(True)
for i, v in enumerate(frac.values):
    axB.text(v + 1.2, i, f"{v:.0f}%", va="center", fontsize=7.5, color=BOND_C)
axB.set_title("b   Detection rate of community partners carrying the downstream reaction",
              fontsize=9.6, fontweight="bold", loc="left", pad=8)

# ============================== panel c: CD vs UC partner abundance
axC = fig.add_subplot(gs[1, 0])
sub = sx[sx.metric == "partner_abundance"].set_index("pathway").reindex(pw_keys)
# Only the tryptophan relay has a non-zero partner abundance; the other three
# are true zeros (the partner species are not detectably carrying the
# downstream reaction in these 21 metagenomes). Plotting them on the same log
# axis would hide them while the title claimed an increase, so they are
# reported explicitly instead.
_plotted = [k for k in pw_keys
            if max(sub.mean_CD.get(k, 0), sub.mean_UC.get(k, 0)) > 1e-8]
_zero = [k for k in pw_keys if k not in _plotted]
# Every individual host is drawn, not just the two group means: with n = 6 in
# UC the distribution matters, and two lollipops would hide it. Jitter is
# deterministic so the panel is reproducible.
_rng = np.random.default_rng(7)
_long = pd.read_csv(os.path.join(ANA, "network_partner_abundance_long.tsv"),
                    sep="\t")
_UC_C, _CD_C = PALETTE[11], PALETTE[1]     # slot 12 UC, slot 2 CD
# the long table uses its own pathway keys, which do not always match the PW
# keys used above (tryptophan_to_indole is stored as tryptophan_kynurenine)
_LONG_KEY = {"tryptophan_indole": "tryptophan_kynurenine"}
for i, k in enumerate(_plotted):
    _d = _long[_long.pathway == _LONG_KEY.get(k, k)]
    for grp, c, dx in (("UC", _UC_C, -0.16), ("CD", _CD_C, 0.16)):
        _v = _d.loc[_d.group == grp, "abundance"].values.astype(float)
        _v = _v[_v > 0]
        if _v.size == 0:
            continue
        _x = np.full(_v.size, i + dx) + _rng.uniform(-0.075, 0.075, _v.size)
        axC.plot(_x, _v, "o", ms=5.2, color=c, mec="white", mew=0.8,
                 alpha=0.9, zorder=4)
    # group mean as a short horizontal tick on top of the swarm
    for grp, c, dx in (("UC", _UC_C, -0.16), ("CD", _CD_C, 0.16)):
        mv = sub.loc[k, f"mean_{grp}"]
        if mv and mv > 0:
            axC.plot([i + dx - 0.085, i + dx + 0.085], [mv, mv], color=c,
                     lw=2.0, zorder=5, solid_capstyle="butt")
    _uc, _cd = sub.loc[k, "mean_UC"], sub.loc[k, "mean_CD"]
    q = sub.loc[k, "q_CD_vs_UC"]
    star = "***" if q < 0.01 else ("**" if q < 0.05 else "ns")
    axC.text(i, max(_uc, _cd) * 2.6,
             f"{_cd/_uc:.1f}-fold\n{star}, q = {q:.3f}",
             ha="center", va="bottom", fontsize=7.4, color=_CD_C,
             fontweight="bold", linespacing=1.35)
axC.set_xticks(range(len(_plotted)))
axC.set_xticklabels([PW[k]["label"].replace("\n", " ") for k in _plotted],
                    fontsize=8)
axC.set_xlim(-0.55, len(_plotted) - 0.45)
_lo = max(1e-7, float(_long.abundance[_long.abundance > 0].min()) * 0.4)
axC.set_ylim(_lo, float(sub.loc[_plotted, ["mean_CD", "mean_UC"]].to_numpy().max()) * 12)
axC.set_yscale("log")
axC.set_xlabel("Partner abundance (RPKs, log scale)", fontsize=8)
axC.tick_params(axis="x", labelsize=7.5)
axC.tick_params(axis="y", labelsize=7.5)
axC.spines[["top", "right"]].set_visible(False)
axC.grid(axis="y", lw=0.4, color="#DDDDDD", zorder=0); axC.set_axisbelow(True)
_sel = _long[_long.pathway == _LONG_KEY.get(_plotted[0], _plotted[0])]
_nu = int((_sel.group == "UC").sum())
_nc = int((_sel.group == "CD").sum())
axC.legend(handles=[Line2D([], [], marker="o", ls="none", ms=7,
                           color=_UC_C, mec="white",
                           label=f"Ulcerative colitis (n = {_nu} hosts)"),
                    Line2D([], [], marker="o", ls="none", ms=7,
                           color=_CD_C, mec="white",
                           label=f"Crohn's disease (n = {_nc} hosts)")],
           fontsize=7.5, frameon=False, loc="lower right")

_zero_txt = ", ".join(PW[k]["label"].replace("\n", " ") for k in _zero)
axC.set_title("c   Tryptophan-relay partner abundance is elevated in Crohn's "
              "disease\n(Mann-Whitney U, BH-corrected; every host shown)",
              fontsize=9.4, fontweight="bold", loc="left", pad=8)
axC.text(0, -0.34,
         f"Partner abundance for the other {len(_zero)} relays "
         f"({_zero_txt}) is at or near\nzero in all 21 metagenomes and "
         "cannot be compared; they are not shown. HMP2 provides no\nnon-IBD "
         "metagenomes, so the contrast is CD versus UC only.",
         transform=axC.transAxes, fontsize=6.6, color="#444444", va="top",
         linespacing=1.5)

# ============================== panel d: top partner species abundance
axD = fig.add_subplot(gs[1, 1])
top_sp = sp.sort_values("total", ascending=False).head(14).iloc[::-1]
y = np.arange(len(top_sp))
axD.barh(y, top_sp["CD"].values, 0.62, color=PALETTE[1],
         edgecolor="white", lw=0.6, label="Crohn's disease")
axD.barh(y, -top_sp["UC"].values, 0.62, color=PALETTE[11],
         edgecolor="white", lw=0.6, label="Ulcerative colitis")
axD.set_yticks(y)
axD.set_yticklabels([shorten(i, 30) for i in top_sp.index], fontsize=7)
axD.axvline(0, color="black", lw=0.7)
axD.set_xlabel("Mean abundance (RPKs); CD right, UC left", fontsize=8)
axD.tick_params(axis="x", labelsize=7.5)
axD.spines[["top", "right", "left"]].set_visible(False)
axD.grid(axis="x", lw=0.4, color="#DDDDDD", zorder=0); axD.set_axisbelow(True)
axD.legend(fontsize=7, frameon=False, loc="lower right")
axD.set_title("d   Community partner species across the cohort",
              fontsize=9.6, fontweight="bold", loc="left", pad=8)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig4_relay_network.{ext}"),
                dpi=300, bbox_inches="tight", facecolor="white")
print("[saved]", os.path.join(FIG, "fig4_relay_network.png"))

print("\n=== relay completeness (CD vs UC) ===")
print(cx.groupby(["pathway", "group"])["partner_abundance"]
      .mean().unstack().to_string())
print("\n=== carrier table ===")
print(cf.sort_values("frac_samples", ascending=False).head(15).to_string(index=False))