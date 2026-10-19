#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GLOBAL figure palette - single source of truth for main AND supplementary
figures.
=============================================================================
The user specified a 13-colour palette and requires that (a) the palette is
used in the given fixed order, and (b) every figure in the paper - main text
and supplementary alike - draws from the same ordered list so the figures read
as one visual system.

    1 #EE7F31 orange  2 #DC2234 red       3 #1F9A49 green
    4 #7F649D purple  5 #0173B2 blue      6 #C7C428 yellow
    7 #79B549 l.green 8 #4FA988 teal      9 #88554F brown
   10 #1C9D49 d.green 11 #0176B2 deepblue 12 #E57661 salmon
   13 #CA79A6 pink

DUPLICATE-HUE WARNING (important, keep in mind when adding series)
--------------------------------------------------------------------
#1F9A49 / #1C9D49 are both green and #0173B2 / #0176B2 are both blue. Entries
10and 11 therefore must never be the only distinguishing feature between two
adjacent series in the same panel. When a figure needs more than 8 categories,
those slots are used only together with a different marker/linestyle/hatch.

FIXED SEMANTIC ASSIGNMENT (colour threading - one entity keeps one colour in
every figure of the paper, main and supplementary alike):
    1 orange  = UPSTREAM, substrate-generating, "probiotic-supplied"
    2 red     = DOWNSTREAM, host-effect-generating, disease (CD)
    3 green   = PARTNER organism (community member)
    4 purple  = Tryptophan relay / free-form accent / statistical highlight
    5 blue    = Bile-acid relay
    6 yellow  = lactate-SCFA relay
    7 l.green = GABA relay
    8 teal    = TMA/TMAO relay
    9 brown   = neutral reference / core-orthogroup control
   10 d.green= combination / rescue / validated (slot 3 + hatch or marker)
   11 d.blue = cross-cohort / secondary series (slot 5 + hatch or marker)
   12 salmon = clinical metabolome, ulcerative colitis
   13 pink   = clinical metabolome, ulcerative-colitis secondary series
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# ----------------------------------------------------------------- the palette
PALETTE = [
    "#EE7F31",  #  1 orange
    "#DC2234",  #  2 red
    "#1F9A49",  #  3 green
    "#7F649D",  #  4 purple
    "#0173B2",  #  5 blue
    "#C7C428",  #  6 yellow
    "#79B549",  #  7 light green
    "#4FA988",  #  8 teal
    "#88554F",  #  9 brown
    "#1C9D49",  # 10 dark green   (near-duplicate of 3)
    "#0176B2",  # 11 deep blue    (near-duplicate of 5)
    "#E57661",  # 12 salmon
    "#CA79A6",  # 13 pink
]

# ------------------------------------------------- fixed semantic roles
UP_C = PALETTE[0]        # 1  upstream / substrate-generating
DOWN_C = PALETTE[1]      # 2  downstream / host-effect-generating
BOND_C = PALETTE[2]      # 3  partner organism / combination
TRP_C = PALETTE[3]       # 4  tryptophan relay
BA_C = PALETTE[4]        # 5  bile-acid relay
LAC_C = PALETTE[5]       # 6  lactate-SCFA relay
GABA_C = PALETTE[6]      # 7  GABA relay
TMA_C = PALETTE[7]       # 8  TMA/TMAO relay
NEUTRAL_C = PALETTE[8]   # 9  reference / core-OG control
RESCUE_C = PALETTE[9]    # 10 combination / rescue  (pairs with slot 3)
SECOND_C = PALETTE[10]   # 11 secondary series     (pairs with slot 5)
UC_C = PALETTE[11]       # 12 ulcerative colitis
CD_C = PALETTE[12]       # 13 Crohn's-disease secondary

NEUTRAL = "#9E9E9E"      # tree branches / grid, deliberately outside the palette
ABSENT = "#F2F2F2"

# --------------------------------------------- relay pathway colour lookup
# Fixed order: lactate-SCFA, bile acid, TMA/TMAO, GABA, tryptophan.
# Slots 6, 5, 8, 7, 4 - taken in palette order, one per pathway, never reused
# by any other entity in the paper.
RELAY_COLORS = {
    "lactate_scfa":  PALETTE[5],
    "bile_acid":     PALETTE[4],
    "tma_tmao":      PALETTE[7],
    "gaba_glutamate": PALETTE[6],
    "tryptophan":    PALETTE[3],
}
RELAY_ORDER = ["lactate_scfa", "bile_acid", "tma_tmao",
               "gaba_glutamate", "tryptophan"]

DISEASE_COLORS = {"CD": PALETTE[1], "UC": PALETTE[11], "nonIBD": PALETTE[8]}


# ----------------------------------------------------------------- helpers
def pick(*indices):
    """Palette colours for 1-based positions (index 1 = first colour)."""
    return [PALETTE[(int(i) - 1) % len(PALETTE)] for i in indices]


def seq(n):
    """The first n palette colours, strictly in the user's fixed order."""
    return PALETTE[:n]


def present_cmap(color=UP_C, name="present"):
    """Discrete two-step single-hue ramp for a binary trait matrix.

    Single hue at two lightness levels instead of a red/green pair, so the
    panel stays readable under deuteranopia (figure-style 4.5 / 8).
    """
    return ListedColormap([ABSENT, color])


def diverging_cmap(name="div"):
    """Diverging ramp for z-scores, anchored on the NEUTRAL grey so that zero
    is visibly neutral. Bounded by the two palette slots that carry the
    downstream/partner meaning (2 and 5) rather than by a stock matplotlib map.
    """
    return LinearSegmentedColormap.from_list(
        name, [PALETTE[4], "#C8DCEA", "#F5F5F5", "#F2C7C9", PALETTE[1]])


def sequential_cmap(color=PALETTE[4], name="seq"):
    """Single-hue sequential ramp for a continuous quantity."""
    return LinearSegmentedColormap.from_list(name, ["#EFF4F8", color])


def discrete(n, name="disc"):
    """n categorical colours in palette order (n <= 8 recommended)."""
    from matplotlib.colors import BoundaryNorm
    cmap = ListedColormap(PALETTE[:n])
    return cmap, BoundaryNorm(range(n + 1), n)


def apply_style(sizes=(9, 8, 7)):
    """Font ladder + outward ticks + frameless look, identical for every figure."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
        "axes.titlesize": sizes[0],
        "axes.labelsize": sizes[0],
        "legend.fontsize": sizes[1],
        "xtick.labelsize": sizes[2],
        "ytick.labelsize": sizes[2],
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


BASE = r"D:\Probiotics_Pangenome_v5"
FIG_MAIN = os.path.join(BASE, "figures", "main_text")
FIG_SUPP = os.path.join(BASE, "figures", "supplementary")
for _d in (FIG_MAIN, FIG_SUPP):
    os.makedirs(_d, exist_ok=True)

if __name__ == "__main__":
    print("palette order:")
    for i, c in enumerate(PALETTE, 1):
        print(f"  {i:2d}  {c}")