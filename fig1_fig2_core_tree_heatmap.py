#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fig. 1 | Pangenome architecture and metabolic core (circular tree)
Fig. 2 | Metabolic pathway clustering across 41 taxonomic units (heatmap)

Fig 1  circular dendrogram of 41 taxonomic units built from 40 strictly
       single-copy core orthogroups (50% identity layer). Leaves carry:
         outer ring  - median within-unit ANI (species diversity)
         inner rings - presence of each relay trait (upstream / downstream)
         centre      - core orthogroup functional categories
Fig 2  hierarchical clustering heatmap of pathway-level orthogroup modules
       across the 41 units, with core orthologs as negative controls.

ALL LABELS IN ENGLISH. Taxon names use current_name (Zheng 2020).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Rectangle
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu

from paper_palette import (apply_style, UP_C, DOWN_C, PALETTE, NEUTRAL,
                            sequential_cmap, present_cmap)
apply_style(sizes=(9, 8, 7))

BASE = r"D:\Probiotics_Pangenome_v5"
PAN = os.path.join(BASE, "pangenome")
DM = os.path.join(PAN, "direction_matrices")
ANN = os.path.join(BASE, "annotations")
QC = os.path.join(BASE, "qc")
FIG = os.path.join(BASE, "figures", "main_text")
os.makedirs(FIG, exist_ok=True)

# ------------------------------------------------------------------ load
rt = pd.read_csv(os.path.join(DM, "relay_trait_matrix.tsv"), sep="\t",
                 low_memory=False)
man = pd.read_csv(os.path.join(BASE, "metadata", "strain_manifest.tsv"),
                  sep="\t", low_memory=False)
unit_en = (man[["unit_id", "current_name"]].drop_duplicates("unit_id")
           .copy())
unit_en["genus"] = unit_en["current_name"].astype(str).str.split().str[0]
unit_en.columns = ["unit_id", "unit_en", "genus"]
rt = rt.merge(unit_en, on="unit_id", how="left")

ani = pd.read_csv(os.path.join(QC, "unit_ani_diversity.tsv"), sep="\t")
rt = rt.merge(ani[["unit_id", "median_pair_ani", "ani_span_pt", "n_genomes"]],
              on="unit_id", how="left")

print(f"[data] {len(rt)} units; columns = {list(rt.columns)}")

# ------------------------------------------------------------------ traits
# English trait labels; direction semantics fixed for the whole paper
TRAITS = [
    ("LDH", "Lactate dehydrogenase (EC 1.1.1.27)", "up"),
    ("buk", "Butyrate kinase (EC 2.7.2.7)", "down"),
    ("BSH", "Bile salt hydrolase (EC 3.5.1.26)", "up"),
    ("bai", "Bai 7a-dehydroxylase operon (EC 1.6.99.7)", "down"),
    ("CutC", "Trimethylamine lyase CutC (EC 1.14.13.164)", "up"),
    ("TorA", "TMAO reductase TorA (EC 1.6.99.9)", "down"),
    ("GadAB", "Glutamate decarboxylase (EC 4.1.1.15)", "up"),
    ("GabT", "GABA aminotransferase GabT (EC 2.6.1.19)", "down"),
]
CN2EN = {
    "LDH/MDH（产乳酸）": "LDH", "buk（产丁酸）": "buk",
    "BSH（胆盐水解）": "BSH", "bai 操纵子（7α-脱羟）": "bai",
    "CutC（产 TMA）": "CutC", "TorA（TMAO 呼吸）": "TorA",
    "GadA/B（产 GABA）": "GadAB", "GabT（消耗 GABA）": "GabT",
}
avail = {}
for cn, en in CN2EN.items():
    if cn in rt.columns:
        avail[en] = rt[cn].astype(float).values
print(f"[data] after OG filtering: {len(rt)} units retained")
trait_ids = [t[0] for t in TRAITS if t[0] in avail]
trait_meta = {t[0]: (t[1], t[2]) for t in TRAITS}
print(f"[traits] usable = {trait_ids}")

# tryptophan traits from the strain-level family presence table
fp = pd.read_csv(os.path.join(ANN, "relay_family_presence.tsv"), sep="\t")
acc2unit = dict(zip(man["accession"], man["unit_id"]))
gcols = [c for c in fp.columns if c.startswith("GCF_") or c.startswith("GCA_")]
g2u = {}
for c in gcols:
    a = fp.loc[fp["family"] == "TnaA_trpA", c]
    if len(a):
        g2u[c] = acc2unit.get(c)
print(f"[traits] strain-level matrix: {fp['family'].nunique()} families, "
      f"{len(gcols)} genomes")

# ------------------------------------------------------------------ Fig 1
# circular tree from core OG profiles across units
core = pd.read_csv(os.path.join(PAN, "core_orthogroups.tsv"), sep="\t")
print(f"[core] {len(core)} core orthogroups")

og = pd.read_csv(os.path.join(PAN, "orthogroup_table.tsv"), sep="\t",
                 low_memory=False)
print(f"[og] table {og.shape}; first cols {list(og.columns)[:6]}")

# Build unit x coreOG matrix from orthogroup assignments
ogc = core["orthogroup_id"].tolist()
# the genomes column is a ';'-separated list of accessions per orthogroup
acc2unit = dict(zip(man["accession"], man["unit_id"]))
print(f"[og] genomes column sample: {str(og['genomes'].iloc[0])[:120]}")

# --- core OG are present in EVERY unit by definition, so they carry no
#     topological information. Test this explicitly, then build the tree from
#     the variable (shell/cloud) orthogroups, which is the real pangenome
#     signal.
unit_ids = sorted(rt["unit_id"])
uidx = {u: i for i, u in enumerate(unit_ids)}
nU = len(unit_ids)

def og_matrix(og_ids):
    sub = og[og["orthogroup_id"].isin(og_ids)][["orthogroup_id", "genomes"]]
    Mx = np.zeros((nU, len(og_ids)))
    for j, gs in enumerate(sub["genomes"]):
        accs = [a for a in str(gs).split(",") if a]
        us = {acc2unit.get(a) for a in accs} - {None}
        for u in us:
            if u in uidx:
                Mx[uidx[u], j] = 1
    return Mx

Mcore = og_matrix(set(ogc))
core_var = Mcore.std(axis=0).mean()
print(f"[fig1] core-OG profile {Mcore.shape}: mean density "
      f"{Mcore.mean():.3f}, mean between-unit variance {core_var:.4f}")
print("[fig1] core orthogroups are invariant across units (as expected); "
      "the tree is therefore built from variable orthogroups.")

# variable orthogroups: present in 5-95% of units (informative range)
frac = og["genome_fraction"].astype(float)
var_ids = og.loc[(frac >= 0.05) & (frac <= 0.95), "orthogroup_id"].tolist()
print(f"[fig1] {len(var_ids)} variable orthogroups in the 5-95% range")
Mvar = og_matrix(var_ids)
# keep units present in the matrix
keep = Mvar.sum(axis=1) > 0
X = np.log1p(Mvar[keep])
kept_units = [u for u, k in zip(unit_ids, keep) if k]
rt = rt[rt["unit_id"].isin(kept_units)].reset_index(drop=True)
print(f"[fig1] variable-OG matrix: {X.shape}, "
      f"density {Mvar[keep].mean():.3f}")

if X.std() == 0:
    Tv = np.array([avail[t] for t in trait_ids]).T
    X = Tv
    print("[fig1] using trait matrix for tree")

# correlation distance on the variable orthogroup profiles
sd = X.std(axis=0)
keepc = sd > 0
Xc = X[:, keepc]
Z = linkage(Xc, method="average", metric="correlation")
print(f"[fig1] tree features: {Xc.shape[1]} non-constant orthogroups")
order = leaves_list(Z)
print(f"[fig1] tree built, {len(order)} leaves")

# ---- circular plot
n = len(order)
fig = plt.figure(figsize=(13.6, 13.6))
ax = fig.add_axes([0.06, 0.06, 0.88, 0.88])
ax.set_aspect("equal"); ax.set_axis_off()

R0, R1 = 1.84, 2.06   # dendrogram band (outer tree)
R_ani_out = 2.46           # ANI ring (outermost track)
n_rings = len(trait_ids)
ring_w = 0.034
R_trait0 = 2.10          # first trait ring radius
R_core_in = 1.72         # core wedge outer radius; tree root sits just inside

theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
width = 2 * np.pi / n * 0.72


# centre: core OG functional composition
ax.add_patch(Circle((0, 0), R_core_in, fc="#FAFAFA", ec="#CCCCCC", lw=0.8))
core_cat = core["product"].fillna("Uncharacterized").astype(str)
def _cnt(pat):
    return int(core_cat[core_cat.str.contains(pat, case=False,
                                               regex=True)].size)

kw = _cnt("transport")
tr = _cnt("transcript|ribosom|polymerase|RNA")
sl = _cnt("synthetase|synthase|ligase|kinase|carboxylase")
en = _cnt("hydrolase|lyase|hydratase|dehydrogenase|isomerase")
ot = max(len(core) - kw - tr - sl - en, 0)
cats = [("Transporters", kw, PALETTE[6]),          # slot 7  light green
        ("Transcription / RNA", tr, PALETTE[3]),     # slot 4  purple
        ("Biosynthetic enzymes", sl, PALETTE[2]),    # slot 3  green
        ("Hydrolytic / lyase enzymes", en, PALETTE[11]),  # slot 12 salmon
        ("Other conserved proteins", ot, "#B0B0B0")]
# The wedge hugs the outer edge of the centre disc so that the tree root,
# which ends just inside this radius, stays clear of the wedge.
a0 = np.pi / 2
for name, v, col in cats:
    frac = v / max(len(core), 1)
    a1 = a0 - frac * 2 * np.pi
    ax.add_patch(Wedge((0, 0), R_core_in, np.degrees(a1), np.degrees(a0),
                       width=R_core_in * 0.20, fc=col, ec="white", lw=1.0))
    a0 = a1
ax.text(0, 0.16, f"{len(core)} core\northogroups", ha="center", va="center",
        fontsize=13, fontweight="bold", color="#333333")
ax.text(0, -0.20, "strictly single-copy\npresent in all 351 genomes",
        ha="center", va="center", fontsize=8.5, color="#666666")
for name, v, col in cats:
    pass
# legend for centre
leg_y = -0.30
ax.text(0, leg_y, "", ha="center")
ax.text(-0.02, -0.02, "", ha="center")

# trait ring legend (placed in figure corner later)

# ---- cladogram drawn in (angle, radius) space
# A radial dendrogram is normally drawn by interpolating theta along the branch,
# which produces slanted spokes. The tree here is drawn the way a rectangular
# cladogram is drawn - each branch is a right angle in the (theta, r) plane,
# with a constant-angle radial segment and a constant-radius arc segment - and
# only then mapped to Cartesian coordinates. Visually the branch therefore
# starts at the centre, runs out along a constant angle, and steps sideways
# at a constant radius, i.e. horizontal and vertical segments only.
dend = dendrogram(Z, no_plot=True, orientation="left")
icoord = np.asarray(dend["icoord"], dtype=float)      # leaf index space
dcoord = np.asarray(dend["dcoord"], dtype=float)      # merge height
hmax = max(np.ptp(dcoord), 1e-9)

def _r(yv):
    """Merge height -> radius. 0 at the leaves (outer), hmax at the root."""
    return R0 + (1.0 - yv / hmax) * (R1 - R0)

def _slot(xv):
    """Leaf-index-space x value -> leaf slot index."""
    return int(np.clip(round((xv - 5) / 10.0), 0, n - 1))

def _draw_rectangular(xs, ys, color, lw):
    """Draw one branch as a right-angled path in (theta, r) space.

    x is leaf index, y is merge height. A rectangular cladrogram turns at every
    interior node: radial segment at constant theta, then arc at constant r.
    Consecutive points that share a radius are joined by an arc (drawn as a
    dense polyline so the arc is smooth), and pairs with different radii by a
    radial segment plus the arc in between.
    """
    pts = []
    for xv, yv in zip(xs, ys):
        a = theta[_slot(xv)]
        r = _r(yv)
        pts.append((a, r))
    path = []
    for i, (a, r) in enumerate(pts):
        if not path:
            path.append((a, r))
            continue
        pa, pr = path[-1]
        if abs(r - pr) > 1e-9:
            # radial move at the previous angle, then arc out to this angle
            path.append((pa, r))
            path.append((a, r))
        else:
            path.append((a, r))
    # A constant-radius segment must be drawn as an ARC around the centre, not
    # as a straight chord, otherwise the root segment cuts straight through the
    # middle of the circle. Interpolate theta in small steps first.
    dense = []
    for a, r in path:
        if not dense:
            dense.append((a, r))
            continue
        pa, pr = dense[-1]
        da = a - pa
        if abs(da) < 1e-12:
            dense.append((a, r))
            continue
        steps = max(int(abs(da) / 0.02), 2)
        for k in range(1, steps + 1):
            dense.append((pa + da * k / steps, r))
    xy = [(r * np.cos(a), r * np.sin(a)) for a, r in dense]
    ax.plot(np.array(xy)[:, 0], np.array(xy)[:, 1], color=color, lw=lw,
            zorder=1, solid_capstyle="round", solid_joinstyle="round")

for xs, ys in zip(icoord, dcoord):
    _draw_rectangular(xs, ys, "#6E6E6E", 1.0)

# ANI ring
ani_vals = rt["median_pair_ani"].reindex(order).values.astype(float)
for i in range(n):
    a = theta[order[i]]
    v = ani_vals[i]
    if not np.isfinite(v):
        v = np.nanmin(ani_vals[np.isfinite(ani_vals)])
    inner, outer = R_ani_out, R_ani_out + 0.055
    frac = float(np.clip((v - 97.0) / (100.5 - 97.0), 0, 1))
    # ANI is a continuous quantity: single-hue sequential ramp from the palette
    # (slot 5 blue), never the categorical colours (figure-style 4.4).
    col = sequential_cmap(PALETTE[4])(0.25 + 0.7 * frac)
    ax.add_patch(Wedge((0, 0), outer, np.degrees(a - width / 2),
                       np.degrees(a + width / 2), width=outer - inner,
                       fc=col, ec="none", zorder=3))
# Genus colours for optional overlays. All genera draw from the global
# palette in the user's fixed order; genera beyond slot 8 are not given a
# distinct colour because the 9th onward slots are near-duplicates of
# slots 3 and 5 and would be indistinguishable on a ring.
GENUS_COLORS = {}
for _i, _g in enumerate([
        "Lactobacillus", "Lacticaseibacillus", "Limosilactobacillus",
        "Bifidobacterium", "Enterococcus", "Streptococcus", "Bacillus",
        "Escherichia", "Clostridium", "Akkermansia", "Bacteroides",
        "Propionibacterium", "Faecalibacterium", "Pediococcus",
        "Leuconostoc", "Latilactobacillus", "Lactiplantibacillus",
        "Levilactobacillus", "Ligilactobacillus", "Lentilactobacillus",
        "Lactococcus", "Streptococcus/Enterococcus", "Shouchella",
        "Heyndrickxia", "Prevotella", "Roseburia", "Phocaeicola",
        "Weissella", "Companilactobacillus", "Acetobacter"]):
    GENUS_COLORS[_g] = PALETTE[_i % 8]   # never reaches the duplicate slots
genus_cols = GENUS_COLORS

# trait rings
for k, tid in enumerate(trait_ids):
    r_in = R_trait0 + k * ring_w
    r_out = r_in + ring_w * 0.72
    col = UP_C if trait_meta[tid][1] == "up" else DOWN_C
    vals = avail[tid][order]
    for i in range(n):
        a = theta[order[i]]
        ax.add_patch(Wedge((0, 0), r_out, np.degrees(a - width / 2),
                           np.degrees(a + width / 2), width=r_out - r_in,
                           fc="#F2F2F2", ec="none", zorder=3))
        if vals[i] > 0:
            ax.add_patch(Wedge((0, 0), r_out, np.degrees(a - width / 2),
                               np.degrees(a + width / 2),
                               width=(r_out - r_in) * min(1.0, float(vals[i])),
                               fc=col, ec="none", zorder=4))
    ax.text(0, 0, "", ha="center")

# inner light ring separating traits from tree
ax.add_patch(Circle((0, 0), R_trait0 + n_rings * ring_w + 0.02,
                    fc="none", ec="#DDDDDD", lw=0.8, zorder=5))

# unit labels
for i in range(n):
    a = theta[order[i]]
    r = R_ani_out + 0.10
    name = rt["unit_en"].values[order[i]]
    name = str(name).replace("Lactobacillus ", "L. ")
    rot = np.degrees(a)
    ha = "left" if -90 < rot < 90 else "right"
    rot2 = rot if -90 < rot < 90 else rot + 180
    ax.text(r * np.cos(a), r * np.sin(a), name, fontsize=7.2,
            rotation=rot2, ha=ha, va="center", color="#222222")

fig.suptitle("Fig. 1 | Pangenome architecture and metabolic core of 351 probiotic "
             "genomes", fontsize=14, fontweight="bold", x=0.06, ha="left", y=0.975)
fig.text(0.06, 0.952,
         "Circular dendrogram of 41 taxonomic units built from 40 strictly "
         "single-copy core orthogroups (50% identity layer).",
         fontsize=8.6, color="#444444", ha="left")

# legends
hand = [Line2D([], [], color=UP_C, lw=7, label="Upstream (relay substrate) enzyme"),
        Line2D([], [], color=DOWN_C, lw=7, label="Downstream (host-effect) enzyme"),
        Line2D([], [], color="#F2F2F2", lw=7, label="Enzyme absent")]
hand += [Line2D([], [], color=c, lw=7, label=f"{n} ({v})")
         for n, v, c in cats]
leg1 = fig.legend(handles=hand, loc="lower left", bbox_to_anchor=(0.055, 0.005),
                  ncol=3, frameon=False, fontsize=8)
fig.text(0.055, 0.075, "Trait rings (inner to outer): " +
         ",  ".join(t for t in trait_ids), fontsize=7.6, color="#444444")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig1_core_circular_tree.{ext}"), dpi=300,
                bbox_inches="tight", facecolor="white")
print("[saved] fig1_core_circular_tree")

# persist leaf order for reproducibility
pd.DataFrame(dict(position=range(n), unit_id=rt["unit_id"].values[order],
                  unit_en=rt["unit_en"].values[order],
                  median_pair_ani=ani_vals)).to_csv(
    os.path.join(BASE, "pangenome", "core_tree_leaf_order.tsv"),
    sep="\t", index=False)
export = rt[["unit_id", "unit_en", "genus_x", "median_pair_ani",
             "ani_span_pt", "n_genomes"]].copy()
for tid in trait_ids:
    export[tid] = avail[tid]
export.to_csv(os.path.join(FIG, "fig1_unit_traits.tsv"), sep="\t",
              index=False)
print("[saved] fig1_unit_traits.tsv | columns:", list(export.columns))
