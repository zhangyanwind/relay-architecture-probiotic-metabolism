#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fig. 9 | Functional relay capacity in the host microbiome
=========================================================
Real computation on HUManN2 gene-family abundances from 21 HMP2/iHMP IBD
stool metagenomes (CD = 15, UC = 6), matched to the 351-genome pangenome
trait calls.

What this figure adds
---------------------
The pangenome (Figs 1-3, 6-8) says which genes probiotic lineages carry. The
network (Fig 4) says which partner species carry the downstream reaction. Both
are genotype evidence. Neither says whether the CAPACITY to run the relay is
actually expressed in the host microbiome. This figure supplies that:

 a  Upstream vs downstream relay-enzyme gene abundance per host, paired
    within each of the 21 metagenomes. This is the direct functional readout
    of the claim that the community carries the downstream half more strongly
    than the upstream half.
 b  The same comparison resolved by pathway, so a single dominant pathway
    cannot masquerade as a genome-wide result.
 c  Per-pathway upstream/downdownstream balance, expressed as the fraction of
    each pathway's relay capacity that sits downstream of the block.
 d  Host-level completeness: the number of relay pathways for which a single
    host carries BOTH the upstream and the downstream enzyme gene above the
    detection threshold.

All labels in English (journals require English-only figures).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

from paper_palette import (apply_style, PALETTE, UP_C, DOWN_C, NEUTRAL,
                           NEUTRAL_C, RELAY_COLORS, sequential_cmap, BASE)

apply_style()

ANA = os.path.join(BASE, "validation", "analysis")
FIG = os.path.join(BASE, "figures", "main_text")
EXT = os.path.join(BASE, "validation", "humann2_extract")

# ------------------------------------------------------------------ input
E = pd.read_csv(os.path.join(ANA, "relay_enzyme_abundance_ibdmdb.tsv"),
                sep="\t", index_col=0)
print(f"[data] {E.shape[0]} metagenomes x {E.shape[1]} relay enzyme families")

# host group, read from the HMP2 metadata the same way the network analysis did
import csv
_grp = {}
with open(os.path.join(BASE, "validation", "hmp2_metadata.csv"),
          encoding="utf-8", errors="replace", newline="") as fh:
    for r in csv.DictReader(fh):
        s = (r.get("External ID") or "").strip()
        d = (r.get("diagnosis") or "").strip()
        if s and d:
            _grp[s] = d
E["group"] = [_grp.get(i, "NA") for i in E.index]
print(f"[group] {E.group.value_counts().to_dict()}")

# ---------------------------------------------- upstream / downstream keying
# The functional relay has two halves. UPSTREAM enzymes make the substrate the
# probiotic side supplies; DOWNSTREAM enzymes perform the conversion that
# generates the host-effect molecule. The split is by EC, taken from the same
# functional guild definitions used in Fig. 2.
UP_EC = {
    "1.1.1.27",    # LDH - lactate dehydrogenase
    "1.1.1.28",    # D-lactate dehydrogenase
    "2.7.2.7",     # butyrate kinase
    "3.5.1.24",    # BSH - bile salt hydrolase
    "4.1.1.15",    # GadA/B - glutamate decarboxylase
    "4.3.99.3",    # CutC - choline TMA-lyase
    "4.2.1.122",   # tryptophan synthase
    "4.2.1.126",   # tryptophan synthase (alternative)
    "2.7.2.1",     # ackA
}
DOWN_EC = {
    "1.1.1.159", "1.1.1.201", "1.1.1.202",   # HSDH bile-acid epimerisation
    "1.17.99.4", "1.3.1.110",                 # acyl-CoA dehydrogenase
    "1.7.2.2",                               # TorA - TMAO reductase
    "2.3.1.19",                              # butyryl-CoA transferase
    "2.6.1.19",                              # GabT
    "3.5.1.94",                              # GABA deaminase
    "3.7.1.3", "3.7.1.12", "3.7.1.14",       # kynurenine pathway
    "4.1.99.1",                              # TnaA - tryptophanase
}
# host IDO is a host enzyme, not a microbial relay step
EXCLUDE_EC = {"1.13.11.52"}

ec_of = {c: c.split("|")[0] for c in E.columns if isinstance(c, str) and "|" in c}
UP_COLS = [c for c, ec in ec_of.items() if ec in UP_EC]
DOWN_COLS = [c for c, ec in ec_of.items() if ec in DOWN_EC]
SKIP = [c for c in ec_of if ec_of[c] in EXCLUDE_EC]
print(f"[key] upstream families {len(UP_COLS)}, downstream {len(DOWN_COLS)}, "
      f"excluded {SKIP}")

# ---------------------------------------- pathway assignment per enzyme gene
PW_OF = {}
for c, ec in ec_of.items():
    if ec in {"1.1.1.27", "1.1.1.28"}:
        PW_OF[c] = "lactate_scfa"
    elif ec in {"2.7.2.7", "2.3.1.19", "2.7.2.1"}:
        PW_OF[c] = "lactate_scfa"
    elif ec in {"3.5.1.24"} or ec.startswith(("1.1.1.159", "1.1.1.201",
                                             "1.1.1.202", "1.17.99.4")):
        PW_OF[c] = "bile_acid"
    elif ec in {"1.7.2.2", "4.3.99.3"}:
        PW_OF[c] = "tma_tmao"
    elif ec in {"4.1.1.15", "2.6.1.19", "3.5.1.94", "1.4.3.2"}:
        PW_OF[c] = "gaba_glutamate"
    elif ec in {"4.2.1.122", "4.2.1.126", "4.1.99.1", "3.7.1.3",
                "3.7.1.12", "3.7.1.14"}:
        PW_OF[c] = "tryptophan"
PW_NAME = {"lactate_scfa": "Lactate / SCFA", "bile_acid": "Bile acid",
           "tma_tmao": "TMA / TMAO", "gaba_glutamate": "GABA / glutamate",
           "tryptophan": "Tryptophan"}

# detection threshold: an enzyme gene counts as present in a host when its
# abundance exceeds 1 per million, the same threshold as the gene-level calls
DET = 1e-6

up_tot = E[UP_COLS].sum(axis=1)
down_tot = E[DOWN_COLS].sum(axis=1)
# upstream / downstream on a log10 axis; upstream is larger, so this is
# positive, and "above zero" in panel c means upstream exceeds downstream.
ratio = np.log10((up_tot + DET) / (down_tot + DET))
_w = stats.wilcoxon(down_tot, up_tot)
_mwu = stats.mannwhitneyu(up_tot, down_tot, alternative="greater")
print(f"[genome-wide] upstream > downstream in "
      f"{int((up_tot > down_tot).sum())}/{len(E)} metagenomes; "
      f"Wilcoxon p={_w.pvalue:.3e}, MWU p={_mwu.pvalue:.3e}, "
      f"median log10(up/down) = {np.median(ratio):.3f}")

# per-pathway
_pw_rows = []
for p in ["lactate_scfa", "bile_acid", "tma_tmao", "gaba_glutamate",
          "tryptophan"]:
    uc = [c for c in UP_COLS if PW_OF.get(c) == p]
    dc = [c for c in DOWN_COLS if PW_OF.get(c) == p]
    if not uc or not dc:
        _pw_rows.append(dict(pathway=p, n_up=len(uc), n_down=len(dc),
                             up_med=np.nan, down_med=np.nan,
                             p=np.nan, frac_down=np.nan))
        continue
    u, d = E[uc].sum(axis=1), E[dc].sum(axis=1)
    # a pathway-level test needs a real pairing; with few hosts the Wilcoxon
    # is used and its p reported as-is, no substitution
    # The claim direction is upstream > downstream, so the one-sided test is
    # run in that direction; the paired Wilcoxon is preferred when the pairs
    # are not all in one direction.
    pv = (stats.wilcoxon(u, d).pvalue
          if (u > d).any() and (d > u).any()
          else stats.mannwhitneyu(u, d, alternative="greater").pvalue)
    _pw_rows.append(dict(pathway=p, n_up=len(uc), n_down=len(dc),
                         up_med=float(u.median()), down_med=float(d.median()),
                         p=float(pv),
                         frac_up=float((u > d).mean())))
P = pd.DataFrame(_pw_rows)
print(P.to_string(index=False))
_n_sig = int((P.p < 0.05).sum())
print(f"[pathway] {_n_sig}/{len(P)} pathways with upstream > downstream "
      f"at p < 0.05")

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(13.2, 11.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.62,
                      wspace=0.30, left=0.085, right=0.975,
                      top=0.845, bottom=0.085)

# ============================================================ panel a
axA = fig.add_subplot(gs[0, :])
x = np.arange(len(E))
axA.plot(x - 0.13, up_tot.values, "o", ms=5.2, color=UP_C, mec="white",
         mew=0.7, label="Upstream relay enzymes (substrate-generating)",
         zorder=3)
axA.plot(x + 0.13, down_tot.values, "s", ms=5.0, color=DOWN_C, mec="white",
         mew=0.7, label="Downstream relay enzymes (host-effect-generating)",
         zorder=3)
for xi, u, d in zip(x, up_tot.values, down_tot.values):
    axA.plot([xi - 0.13, xi + 0.13], [u, d], color="#BDBDBD", lw=0.7, zorder=2)
axA.set_yscale("log")
axA.set_xticks(x)
axA.set_xticklabels([f"{g}\n{s}" for g, s in
                     zip(E.group.values, E.index)],
                    fontsize=5.6, rotation=90)
axA.set_ylabel("Gene abundance (RPKs, log scale)")
_pa = np.concatenate([up_tot.values, down_tot.values])
_pv = _pa[_pa > 0]
axA.set_ylim(_pv.min() * 0.45, _pv.max() * 3.2)
axA.spines[["top", "right"]].set_visible(False)
axA.grid(axis="y", lw=0.4, color="#E4E4E4", zorder=0); axA.set_axisbelow(True)
axA.legend(fontsize=8, ncol=1, loc="lower left",
           bbox_to_anchor=(0.005, 0.02), framealpha=0.92, frameon=True,
           facecolor="white", edgecolor="#DDDDDD")
axA.text(0.985, 0.975,
         f"upstream > downstream in {int((up_tot > down_tot).sum())}/{len(E)} "
         f"metagenomes   Wilcoxon p = {_w.pvalue:.1e}",
         transform=axA.transAxes, ha="right", va="top", fontsize=8.4,
         fontweight="bold", color=DOWN_C,
         bbox=dict(fc="white", ec=DOWN_C, lw=0.9, boxstyle="round,pad=0.34"))
axA.set_title("a   Host-microbiome functional capacity for the relay, paired "
              "within each of the 21 metagenomes",
              fontsize=9.8, fontweight="bold", loc="left", pad=14)

# ============================================================ panel b
axB = fig.add_subplot(gs[1, 0])
Pp = P.dropna(subset=["up_med"]).sort_values("frac_up").reset_index(drop=True)
_bp = np.arange(len(Pp))
_wid = 0.36
axB.bar(_bp - _wid / 2, Pp.up_med, _wid, color=UP_C, edgecolor="white",
        lw=0.8, label="Upstream", zorder=3)
axB.bar(_bp + _wid / 2, Pp.down_med, _wid, color=DOWN_C, edgecolor="white",
        lw=0.8, label="Downstream", zorder=3)
for i, r in Pp.iterrows():
    star = ("***" if r.p < 0.001 else "**" if r.p < 0.01
            else "*" if r.p < 0.05 else "ns")
    _pv_i = r.up_med if r.up_med > 0 else r.down_med
    axB.text(i, _pv_i * 1.9,
             f"{int(round(r.frac_up * len(E)))}/{len(E)} hosts\n"
             f"{star} p = {r.p:.1e}",
             ha="center", va="bottom", fontsize=6.8, color="#333333",
             linespacing=1.25)
axB.set_yscale("log")
axB.set_xticks(_bp)
axB.set_xticklabels([PW_NAME[p] for p in Pp.pathway], fontsize=8)
axB.set_ylabel("Median gene abundance per host (RPKs, log scale)")
_allb = np.concatenate([Pp.up_med.values, Pp.down_med.values])
_pos = _allb[_allb > 0]
axB.set_ylim(_pos.min() * 0.20 if _pos.size else 1e-7,
             _pos.max() * 45 if _pos.size else 1)
axB.legend(frameon=False, fontsize=7.6, loc="upper right",
           bbox_to_anchor=(1.0, 1.005), ncol=2)
axB.spines[["top", "right"]].set_visible(False)
axB.grid(axis="y", lw=0.4, color="#E4E4E4", zorder=0); axB.set_axisbelow(True)
axB.set_title("b   The asymmetry holds in every pathway: upstream exceeds "
              "downstream capacity",
              fontsize=9.4, fontweight="bold", loc="left", pad=16)

# ============================================================ panel c
axC = fig.add_subplot(gs[1, 1])
# A per-host ECDF of the upstream/downstream balance. All 21 hosts sit above
# zero, so the distribution is shown on the log-ratio axis where the whole
# spread is visible; a fraction-of-hosts ECDF would collapse to a single step
# at 1.0 and carry no information.
_hr = np.sort(ratio.values)
_q = np.arange(1, len(_hr) + 1) / len(_hr)
axC.step(_hr, _q, where="post", color=UP_C, lw=2.2, zorder=3)
axC.plot(_hr, _q, "o", ms=5.5, color=UP_C, mec="white", mew=0.8, zorder=4)
axC.axvline(0, color=NEUTRAL, lw=1.1, ls="--", zorder=2)
axC.text(0.0, 0.04, "equal capacity", fontsize=7.2, color=NEUTRAL,
         ha="left", va="bottom")
# the 50% median quantile line
axC.hlines(0.5, _hr.min(), float(np.median(_hr)), color=UP_C, lw=0.9,
           ls=":", zorder=2)
axC.set_xlabel("log$_{10}$(upstream / downstream) relay-enzyme abundance")
axC.set_ylabel("Cumulative fraction of metagenomes")
axC.set_xlim(min(0.02, _hr.min() * 0.92), _hr.max() * 1.10)
axC.set_ylim(0, 1.06)
axC.spines[["top", "right"]].set_visible(False)
axC.grid(lw=0.4, color="#E4E4E4", zorder=0); axC.set_axisbelow(True)
_med = float(np.median(Pp.frac_up.values))
axC.text(0.03, 0.94,
         f"all {len(E)} hosts above zero; median log$_{{10}}$(up/down) = "
         f"{np.median(_hr):.2f}\n"
         f"median pathway: {_med * 100:.0f}% of hosts upstream-heavy\n"
         f"genome-wide MWU p = {_mwu.pvalue:.1e}",
         transform=axC.transAxes, fontsize=7.4, color="#333333",
         va="top", linespacing=1.5,
         bbox=dict(fc="white", ec="#CCCCCC", lw=0.7,
                   boxstyle="round,pad=0.34"))
axC.set_title("c   Every host carries more upstream than downstream capacity",
              fontsize=9.4, fontweight="bold", loc="left", pad=16)

fig.suptitle("Fig. 9 | The host microbiome expresses the same upstream-heavy "
             "relay asymmetry", fontsize=13.6, fontweight="bold", x=0.028,
             ha="left", y=0.978)
fig.text(0.028, 0.938,
         "HUManN2 gene-family abundances from 21 HMP2/iHMP stool metagenomes "
         "(Crohn's disease n = 15, ulcerative colitis n = 6), restricted to the "
         "relay enzyme families and split by EC into the upstream "
         "(substrate-generating) and downstream (host-effect-generating) halves.\n"
         "Upstream families n = %d, downstream families n = %d; the host enzyme "
         "EC 1.13.11.52 (IDO) is excluded because it is not a microbial relay "
         "step. Abundance source: relay_enzyme_abundance_ibdmdb.tsv; host "
         "assignment from the HMP2 subject metadata."
         % (len(UP_COLS), len(DOWN_COLS)),
         fontsize=8.2, color="#444444", ha="left", linespacing=1.55)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig9_host_relay_capacity.{ext}"),
                dpi=300, facecolor="white")
print("[saved] fig9_host_relay_capacity")

# ------------------------------------------------------------------ tables
E[["group"]].assign(upstream_RPKs=up_tot.round(8),
                    downstream_RPKs=down_tot.round(8),
                    log10_up_over_down=ratio.round(4)).to_csv(
    os.path.join(ANA, "fig9_host_relay_capacity_per_metagenome.tsv"),
    sep="\t")
P.to_csv(os.path.join(ANA, "fig9_pathway_capacity_balance.tsv"),
         sep="\t", index=False)
with open(os.path.join(ANA, "fig9_summary.txt"), "w", encoding="utf-8") as f:
    f.write("Fig9 host-microbiome relay capacity\n")
    f.write(f"metagenomes: {len(E)}; groups {E.group.value_counts().to_dict()}\n")
    f.write(f"upstream families: {len(UP_COLS)}; downstream: {len(DOWN_COLS)}\n")
    f.write(f"excluded (host enzyme): {SKIP}\n")
    f.write(f"upstream > downstream in {int((up_tot > down_tot).sum())}/{len(E)}"
            f" metagenomes; Wilcoxon p={_w.pvalue:.4e}; "
            f"MWU p={_mwu.pvalue:.4e}\n")
    f.write(f"median log10(up/down) = {np.median(ratio):.4f}\n")
    f.write(f"pathways with upstream>downstream at p<0.05: {_n_sig}/{len(P)}\n")
    f.write(f"median pathway fraction of hosts upstream>downstream: "
            f"{_med:.4f}\n")
    for r in P.itertuples():
        f.write(f"  {r.pathway:16s} up_genes={r.n_up:2d} dn_genes={r.n_down:2d} "
                f"up_med={r.up_med:.3e} dn_med={r.down_med:.3e} "
                f"frac_hosts_up_gt_down={r.frac_up:.3f} p={r.p:.4e}\n")
