#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fig. 2 | Metabolic pathway modules across 41 probiotic taxonomic units.

Hierarchical clustering heat map built from 24 functional guild modules
covering five metabolic pathways, plus strictly single-copy core
orthogroups as negative controls. Units are clustered by their variable
orthogroup profile, so column order reflects homology rather than taxonomy.

ALL LABELS IN ENGLISH.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list

from paper_palette import (apply_style, PALETTE, UP_C, DOWN_C, NEUTRAL,
                           NEUTRAL_C, RELAY_COLORS, RELAY_ORDER, diverging_cmap)
apply_style(sizes=(9, 8, 7))

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "axes.linewidth": 0.8,
})

BASE = r"D:\Probiotics_Pangenome_v5"
PAN = os.path.join(BASE, "pangenome")
DM = os.path.join(PAN, "direction_matrices")
FIG = os.path.join(BASE, "figures", "main_text")
os.makedirs(FIG, exist_ok=True)

GUILD = {
    "乳酸脱氢酶 LDH/MDH 家族": ("LDH / MDH", "Lactate / SCFA", "up"),
    "D-乳酸脱氢酶": ("D-lactate dehydrogenase", "Lactate / SCFA", "up"),
    "乳酸消旋酶 LarA": ("Lactate racemase LarA", "Lactate / SCFA", "up"),
    "乳酸转运蛋白": ("Lactate transporter", "Lactate / SCFA", "up"),
    "乙酸生成 pta-ackA": ("Acetate synthesis pta-ackA", "Lactate / SCFA", "up"),
    "乙醇生成 adhE": ("Ethanol formation adhE", "Lactate / SCFA", "up"),
    "甲酸生成 pflB": ("Formate formation pflB", "Lactate / SCFA", "up"),
    "甲酸脱氢 fdh/fdn": ("Formate dehydrogenase fdh/fdn", "Lactate / SCFA", "up"),
    "丁酸通路 buk": ("Butyrate kinase (EC 2.7.2.7)", "Lactate / SCFA", "down"),
    "丁酸通路 β-氧化": ("Butyryl-CoA beta-oxidation", "Lactate / SCFA", "down"),
    "丙酸通路 prpD": ("Propionate pathway prpD", "Lactate / SCFA", "down"),
    "丙酸通路 pccB/mmgA": ("Propionate pathway pccB/mmgA", "Lactate / SCFA", "down"),
    "胆盐水解酶 BSH": ("Bile salt hydrolase (EC 3.5.1.26)", "Bile acid", "up"),
    "bai 操纵子（7α-脱羟）": ("Bai 7a-dehydroxylase operon", "Bile acid", "down"),
    "7α-HSDH": ("7a-hydroxysteroid dehydrogenase", "Bile acid", "down"),
    "7β-HSDH": ("7b-hydroxysteroid dehydrogenase", "Bile acid", "down"),
    "12α-HSDH": ("12a-hydroxysteroid dehydrogenase", "Bile acid", "down"),
    "色氨酸合成 trpA/trpB": ("Tryptophan synthase TrpA/TrpB", "Tryptophan", "up"),
    "色氨酸转运 mtr/tnaB": ("Tryptophan transporter mtr/tnaB", "Tryptophan", "up"),
    "色氨酸酶 TnaA（产吲哚）": ("Tryptophanase TnaA (EC 4.1.99.03)", "Tryptophan", "down"),
    "犬尿氨酸酶 KynU": ("Kynureninase KynU", "Tryptophan", "down"),
    "犬尿氨酸通路 KynA/B": ("Kynurenine pathway KynA/KynB", "Tryptophan", "down"),
    "芳香族脱羧酶 TyrDC/Tdc": ("Aromatic amino acid decarboxylase", "Tryptophan", "down"),
    "谷氨酸脱羧酶 Gad": ("Glutamate decarboxylase (EC 4.1.1.15)", "GABA", "up"),
}
PATH_ORDER = ["Lactate / SCFA", "Bile acid", "GABA", "Tryptophan",
              "Core-OG copy-number control"]
# Pathway colours are NOT free choices: they come from the global palette via
# RELAY_COLORS so that a pathway keeps the same colour in every figure of the
# paper (colour threading). Lactate/SCFA = slot 6, bile acid = slot 5,
# GABA = slot 7, tryptophan = slot 4, core-OG control = slot 9.
PC = {"Lactate / SCFA": RELAY_COLORS["lactate_scfa"],
      "Bile acid": RELAY_COLORS["bile_acid"],
      "GABA": RELAY_COLORS["gaba_glutamate"],
      "Tryptophan": RELAY_COLORS["tryptophan"],
      "Core-OG copy-number control": NEUTRAL_C}

guild = pd.read_csv(os.path.join(DM, "unit_functional_guilds.tsv"), sep="\t",
                    low_memory=False)
man = pd.read_csv(os.path.join(BASE, "metadata", "strain_manifest.tsv"),
                  sep="\t", low_memory=False)
u = man[["unit_id", "current_name"]].drop_duplicates("unit_id").copy()
u.columns = ["unit_id", "unit_en"]
guild = guild.merge(u, on="unit_id", how="left").sort_values("unit_id").reset_index(drop=True)
unit_ids = guild["unit_id"].tolist()
unit_en = guild["unit_en"].tolist()
N = len(guild)
print(f"[data] {N} units")

rows, rowlabels, rowgroups = [], [], []
for cn, (en, pw, side) in GUILD.items():
    if cn not in guild.columns:
        continue
    v = pd.to_numeric(guild[cn], errors="coerce").fillna(0).values.astype(float)
    rows.append(np.log1p(v))
    rowlabels.append(f"{en}  ({side})")
    rowgroups.append(pw)
print(f"[guild] {len(rows)} guild modules used")

og = pd.read_csv(os.path.join(PAN, "orthogroup_table.tsv"), sep="\t",
                 low_memory=False)
core = pd.read_csv(os.path.join(PAN, "core_orthogroups.tsv"), sep="\t")
acc2unit = dict(zip(man["accession"], man["unit_id"]))
pos = {u_: i for i, u_ in enumerate(unit_ids)}

sub = og[og["orthogroup_id"].isin(set(core["orthogroup_id"]))]
Mc = np.zeros((N, len(sub)))
for j, gs in enumerate(sub["genomes"]):
    for a in str(gs).split(","):
        uu = acc2unit.get(a)
        if uu in pos:
            Mc[pos[uu], j] = 1
print(f"[core] control matrix {Mc.shape}, density {Mc.mean():.3f}")

rng = np.random.default_rng(42)
core_rows = [Mc.sum(axis=1)]
core_lab = ["Core-OG total copy number"]
sel = rng.choice(Mc.shape[1], 10, replace=False)
for i in sel:
    core_rows.append(Mc[:, i])
    core_lab.append(f"Core OG {int(core['orthogroup_id'].values[i])}")
for lab, v in zip(core_lab, core_rows):
    rows.append(np.asarray(v, float))
    rowlabels.append(lab)
    rowgroups.append("Core-OG copy-number control")

A = np.column_stack(rows)
# z-score each column, but a column with no variance (a core orthogroup is
# present in every genome by definition, so its binary row is constant) must
# collapse to exactly 0. Without this guard the 1e-9 denominator amplifies
# floating-point noise into a +-2 "signal" and the negative control appears
# to vary across units, which would be a fabricated result.
_sd = A.std(axis=0)
_const = _sd < 1e-8
A = np.where(_const, 0.0, (A - A.mean(axis=0)) / np.where(_const, 1.0, _sd))
print(f"[matrix] {A.shape}; {int(_const.sum())} invariant columns "
      f"(negative controls) set to exactly 0")

og["frac"] = og["genome_fraction"].astype(float)
var_ids = og.loc[(og["frac"] >= 0.05) & (og["frac"] <= 0.95),
                "orthogroup_id"].tolist()
subv = og[og["orthogroup_id"].isin(set(var_ids))]
Mv = np.zeros((N, len(subv)))
for j, gs in enumerate(subv["genomes"]):
    for a in str(gs).split(","):
        uu = acc2unit.get(a)
        if uu in pos:
            Mv[pos[uu], j] = 1
sd = Mv.std(axis=0)
Mv = np.log1p(Mv[:, sd > 0])
print(f"[cluster] unit x variable OG {Mv.shape}")
Zu = linkage(Mv, method="average", metric="correlation")
uorder = leaves_list(Zu)

Zrow = linkage(A, method="average", metric="euclidean")
rorder_all = set(leaves_list(Zrow))
rorder = []
for pw in PATH_ORDER:
    idx = {i for i in range(len(rowgroups)) if rowgroups[i] == pw}
    rorder += [i for i in leaves_list(Zrow) if i in idx]
    print(f"[rows] {pw}: {sum(1 for i in rorder if rowgroups[i]==pw)}")

Ad = A[np.ix_(uorder, rorder)]
glabels = [rowgroups[i] for i in rorder]

gb, prev, start = [], None, 0
for i, g in enumerate(glabels):
    if g != prev:
        gb.append((start, i, prev)); prev = g; start = i
gb.append((start, len(glabels), prev))

fig = plt.figure(figsize=(16.0, 12.4))
ax = fig.add_axes([0.275, 0.175, 0.55, 0.685])
# z-score is a signed continuous quantity -> neutral-grey-anchored diverging
# ramp from the global palette, never the stock RdBu (figure-style 4.4).
im = ax.imshow(Ad, aspect="auto", cmap=diverging_cmap("relay_div"),
               vmin=-2.0, vmax=2.0, interpolation="nearest")

# Left pathway strip. The block is drawn INSIDE the axes at negative x with an
# explicit x-limit, so every block has exactly the same width by construction
# (drawing it outside with clip_on=False made the width depend on where the
# row labels happened to end). One thin white separator keeps adjacent
# pathways from bleeding into each other.
BLOCK_W = 6.0
for a, b, g in gb:
    ax.add_patch(Rectangle((-BLOCK_W, a - 0.5), BLOCK_W - 0.6, b - a,
                           fc=PC.get(g, "#999999"), ec="none"))
# separator lines between pathway blocks (full heatmap width)
for a, b, g in gb[1:]:
    ax.axhline(b - 0.5, color="black", lw=1.2)
ax.set_xlim(-BLOCK_W, N - 0.5)
ax.set_ylim(len(rorder) - 0.5, -0.5)

ABBR = [("Lactobacillus ", "L. "), ("Lacticaseibacillus ", "Lc. "),
        ("Limosilactobacillus ", "Lm. "), ("Levilactobacillus ", "Lv. "),
        ("Lactiplantibacillus ", "Lp. "), ("Companilactobacillus ", "Co. "),
        ("Latilactobacillus ", "Lt. ")]
def short_name(x):
    x = str(x)
    for a, b in ABBR:
        x = x.replace(a, b)
    return x

ax.set_xticks(range(N))
ax.set_xticklabels([short_name(unit_en[i]) for i in uorder],
                   rotation=90, fontsize=8)
ax.set_yticks(range(len(rorder)))
ax.set_yticklabels([rowlabels[i] for i in rorder], fontsize=7.6)
ax.set_xlabel("Taxonomic unit, clustered by variable orthogroup profile",
              fontsize=9.5)
ax.tick_params(length=0)
for sp_ in ax.spines.values():
    sp_.set_visible(False)
ax.set_xticks(np.arange(-0.5, N, 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(rorder), 1), minor=True)
ax.grid(which="minor", color="white", lw=0.14)
ax.tick_params(which="minor", length=0)

cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.010)
cb.set_label("Standardised module abundance (z-score)", fontsize=8.5)
cb.ax.tick_params(labelsize=7.5)
cb.outline.set_visible(False)

# Pathway group bar along the BOTTOM edge. The left strip is easy to lose when
# the figure is cropped or scaled, so the same five pathway colours are repeated
# here as a horizontal key in the same fixed order, directly under the columns
# they annotate. Same colours as the left strip - this is a legend, not new data.
axb = fig.add_axes([0.275, 0.098, 0.55, 0.016])
axb.set_xlim(0, N - 0.5)
axb.set_ylim(0, 1)
axb.axis("off")
gw = N / len(PATH_ORDER)
SHORT = {"Lactate / SCFA": "Lactate / SCFA", "Bile acid": "Bile acid",
         "GABA": "GABA", "Tryptophan": "Tryptophan",
         "Core-OG copy-number control": "Core-OG control"}
for k, g in enumerate(PATH_ORDER):
    x0 = k * gw
    axb.add_patch(Rectangle((x0 + 0.12, 0.10), gw - 0.24, 0.80,
                            fc=PC[g], ec="none"))
    axb.text(x0 + gw / 2, 0.5, SHORT[g], ha="center", va="center",
             fontsize=8.4, color="white", fontweight="bold")
axb.text(0, 1.30, "Pathway group", ha="left", va="bottom",
         fontsize=8.5, color="#444444")

fig.suptitle("Fig. 2 | Metabolic pathway modules across 41 probiotic "
             "taxonomic units", fontsize=14.5, fontweight="bold",
             x=0.035, ha="left", y=0.978)
fig.text(0.035, 0.953,
         "Functional guild modules spanning five metabolic pathways, z-scored "
         "per module across units. Core orthogroups provide the copy-number "
         "reference: each is present in every genome by definition, so the "
         "variation in the grey block is copy number, not trait carriage.",
         fontsize=8.5, color="#444444", ha="left")

# The pathway key is now the bar under the heatmap; the remaining free legend
# slot explains the (up)/(down) suffix on every module label, which is the only
# other visual encoding in the panel.
hand = [Line2D([], [], color=UP_C, lw=7, label="upstream / substrate-generating (up)"),
        Line2D([], [], color=DOWN_C, lw=7, label="downstream / host-effect-generating (down)")]
fig.legend(handles=hand, loc="lower left", bbox_to_anchor=(0.035, 0.010),
           ncol=2, frameon=False, fontsize=8.5)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig2_pathway_cluster_heatmap.{ext}"),
                dpi=300, bbox_inches="tight", facecolor="white")
print("[saved] fig2_pathway_cluster_heatmap")

pd.DataFrame(Ad, index=[short_name(unit_en[i]) for i in uorder],
             columns=[rowlabels[i] for i in rorder]).to_csv(
    os.path.join(FIG, "fig2_clustered_matrix.tsv"), sep="\t")
pd.DataFrame(dict(unit_id=[unit_ids[i] for i in uorder],
                  unit_en=[unit_en[i] for i in uorder])).to_csv(
    os.path.join(FIG, "fig2_unit_cluster_order.tsv"), sep="\t", index=False)
print("[saved] tables")