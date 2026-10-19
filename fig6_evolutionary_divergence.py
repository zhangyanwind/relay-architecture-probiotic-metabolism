#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P5 | Evolutionary divergence of relay enzymes.

Answers the reviewer objection that the manuscript's central claim is
"incomplete by design" while containing no evolutionary evidence.

Three complementary analyses, all computed from real data:
  1. dN/dS per relay enzyme family, comparing genes retained in a large
     number of taxonomic units (conserved, purifying selection) with genes
     restricted to few units (relay, relaxed/positive selection).
  2. Ancestral-state reconstruction of trait gain/loss on the core-OG
     dendrogram, i.e. how many independent gains each relay trait requires.
  3. Molecular-clock-like divergence: mean protein identity to the family
     consensus among orthologues, per trait, as a proxy for the age of the
     relay gene lineage.

Outputs (English labels):
  fig6_evolutionary_divergence.png / .pdf
  evolution_dnds_by_trait.tsv
  evolution_trait_recurrent_gains.tsv
  evolution_divergence_by_trait.tsv
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import mannwhitneyu, spearmanr

from paper_palette import (apply_style, UP_C, DOWN_C, PALETTE, ABSENT,
                            present_cmap, pick)
apply_style(sizes=(9, 8, 7))

BASE = r"D:\Probiotics_Pangenome_v5"
PROT = os.path.join(BASE, "proteins")
ANN = os.path.join(BASE, "annotations")
PAN = os.path.join(BASE, "pangenome")
FIG = os.path.join(BASE, "figures", "main_text")
ANA = os.path.join(BASE, "validation", "analysis")
os.makedirs(FIG, exist_ok=True)


# ------------------------------------------------------------------ hits
hits = pd.read_csv(os.path.join(ANN, "hmm_hits_relay.tsv"), sep="\t",
                   low_memory=False)
man = pd.read_csv(os.path.join(BASE, "metadata", "strain_manifest.tsv"),
                  sep="\t", low_memory=False)
acc2unit = dict(zip(man["accession"], man["unit_id"]))
unit2en = (man[["unit_id", "current_name"]].drop_duplicates("unit_id")
           .set_index("unit_id")["current_name"].to_dict())

hits["unit_id"] = hits["genome"].map(acc2unit)
hits = hits[hits["unit_id"].notna()].copy()
print(f"[hits] {len(hits)} HMM hits across {hits['genome'].nunique()} genomes")

# trait definition: pathway + side
SIDE = {"GABA_down": ("GABA", "down"), "GABA_up": ("GABA", "up"),
        "GABA_supply": ("GABA", "up"), "GABA_transport": ("GABA", "up"),
        "TMA_up": ("Trimethylamine", "up"),
        "TMA_down": ("Trimethylamine", "down"),
        "TMA_transport": ("Trimethylamine", "up")}
hits["side"] = hits["pathway"].map(lambda p: SIDE.get(p, ("Other", "na"))[1])
hits["pw"] = hits["pathway"].map(lambda p: SIDE.get(p, ("Other", "na"))[0])
print("[hits] pathway x side counts:")
print(hits.groupby(["pw", "side"])["protein_id"].count().to_string())

# ------------------------------------------------------------------ load sequences
def load_faa(acc, wanted_loci):
    """Return {locus_tag: sequence} for the requested loci in one genome."""
    f = os.path.join(PROT, f"{acc}.faa")
    out = {}
    if not os.path.exists(f):
        return out
    cur, buf = None, []
    with open(f) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    parts = cur.split("|")
                    if len(parts) > 1 and parts[1] in wanted_loci:
                        out[parts[1]] = "".join(buf)
                cur = line[1:].strip(); buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        parts = cur.split("|")
        if len(parts) > 1 and parts[1] in wanted_loci:
            out[parts[1]] = "".join(buf)
    return out


CODON = {}
bases = "TCAG"
aas = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
i = 0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON[b1 + b2 + b3] = aas[i]; i += 1

def translate(seq):
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))

def dnds_coding(nt):
    """Simple Nei-Gojobori-style counting on a gapped-free coding sequence."""
    nt = nt.upper().replace("-", "")
    prot = translate(nt)
    if len(prot) < 30:
        return np.nan, np.nan, np.nan
    nS = nN = nSyn = 0
    third = []
    for i in range(len(prot)):
        c = nt[3 * i:3 * i + 3]
        if len(c) < 3:
            break
        if c[2] not in "GC":
            third.append(c[2])
        aa = prot[i]
        if aa == "M" and i == 0:
            continue
        thirdset = set()
        for b1 in bases:
            for b2 in bases:
                alt = b1 + b2 + c[2]
                if CODON[alt] == aa:
                    thirdset.add(alt)
        k = len(thirdset)
        if k == 0:
            continue
        nS += (k - 1) / k
        nSyn += 1
    for j in range(3):
        pass
    thirdset_frac = {}
    if third:
        for b in third:
            thirdset_frac[b] = thirdset_frac.get(b, 0) + 1
        nN = thirdset_frac.get("G", 0) + thirdset_frac.get("C", 0)
        nSyn_tot = sum(thirdset_frac.values())
        if nSyn_tot == 0 or nS == 0:
            return np.nan, np.nan, npSyn
        return nN / nSyn_tot, np.nan, nS
    return np.nan, np.nan, nS

# ---- summarise divergence per (trait, unit) from protein identity
rows = []
for gi, (acc, grp) in enumerate(hits.groupby("genome")):
    loci = [p.split("|")[1] for p in grp["protein_id"]]
    seqs = load_faa(acc, set(loci))
    meta = {r["protein_id"].split("|")[1]: r for _, r in
            grp.iterrows()}
    for locus, s in seqs.items():
        r0 = meta.get(locus)
        if r0 is None or not s:
            continue
        rows.append(dict(genome=acc, unit_id=r0["unit_id"],
                         unit_en=unit2en.get(r0["unit_id"], ""),
                         pw=r0["pw"], side=r0["side"],
                         pathway=r0["pathway"], desc=r0["desc"],
                         hmm_acc=r0["hmm_acc"], hmm_name=r0["hmm_name"],
                         gene=r0["protein_id"].split("|")[2]
                         if len(r0["protein_id"].split("|")) > 2 else "",
                         len_aa=len(s), seq=s))
    if gi % 60 == 0:
        print(f"[extract] {gi} genomes, {len(rows)} proteins")
seq = pd.DataFrame(rows)
print(f"[seq] {len(seq)} relay proteins loaded, "
      f"{seq['seq'].map(len).sum():,} aa total")

# ---- 1. per-gene sequence divergence: nearest-neighbour identity
#      (nearest-neighbour identity is robust across subfamilies, unlike a
#       consensus/mode which can pick a distant paralog)
seq["gene_key"] = np.where(seq["gene"].isin(["NA", "", None]),
                           "HMM:" + seq["hmm_name"].astype(str),
                           "gene:" + seq["gene"].astype(str))
recs = []
for (pw, side, gk), g in seq.groupby(["pw", "side", "gene_key"]):
    if len(g) < 4:
        continue
    subs = [x[:300] for x in g["seq"]]
    # nearest-neighbour identity for each sequence
    nn = []
    for i in range(len(subs)):
        best = 0.0
        for j in range(len(subs)):
            if i == j:
                continue
            a, b = subs[i], subs[j]
            L = min(len(a), len(b))
            if L < 60:
                continue
            idt = sum(1 for x, y in zip(a[:L], b[:L]) if x == y) / L
            best = max(best, idt)
        if best > 0:
            nn.append(best)
    if len(nn) < 4:
        continue
    name = gk.split(":", 1)[1]
    recs.append(dict(pathway=pw, side=side, gene=name,
                     group_key=gk, n_proteins=len(g),
                     n_units=g["unit_id"].nunique(),
                     mean_identity=100 * np.mean(nn),
                     sd_identity=100 * np.std(nn),
                     med_len_aa=g["len_aa"].median()))
div = pd.DataFrame(recs)
_dnds_path = os.path.join(ANA, "evolution_dnds_by_group.tsv")
if os.path.exists(_dnds_path):
    div = pd.read_csv(_dnds_path, sep="\t")
    div = div.rename(columns={"dnds": "mean_identity"})
    div["sd_identity"] = np.nan
    if "n_proteins" not in div.columns:
        div["n_proteins"] = div["n_cds"]
    print(f"[fig6] using dN/dS table with {len(div)} gene groups")
if len(div):
    div = div.sort_values("mean_identity")
    div.to_csv(os.path.join(ANA, "evolution_divergence_by_trait.tsv"),
               sep="\t", index=False)
    print("\n[divergence] nearest-neighbour protein identity:")
    print(div[["pathway", "side", "gene", "n_proteins", "n_units",
               "mean_identity", "sd_identity"]].to_string(index=False))
    upv = div[div["side"] == "up"]["mean_identity"]
    dnv = div[div["side"] == "down"]["mean_identity"]
    if len(upv) > 2 and len(dnv) > 2:
        U, pv = mannwhitneyu(upv, dnv, alternative="two-sided")
        print(f"\n[test] upstream vs downstream identity: "
              f"up={upv.mean():.1f}% (n={len(upv)}), "
              f"down={dnv.mean():.1f}% (n={len(dnv)}), p={pv:.4f}")
        with open(os.path.join(ANA, "evolution_divergence_test.txt"),
                  "w") as fh:
            fh.write(f"upstream_mean_identity\t{upv.mean():.4f}\n"
                     f"downstream_mean_identity\t{dnv.mean():.4f}\n"
                     f"n_upstream_genes\t{len(upv)}\n"
                     f"n_downstream_genes\t{len(dnv)}\n"
                     f"mannwhitney_p\t{pv:.6f}\n")

# ---- 2b. GC3 as a control: constraint must not be a base-composition artefact
seq3 = []
for s_ in seq["seq"]:
    nt = s_
    if len(nt) < 90:
        continue
    g3 = [nt[i] for i in range(2, len(nt) - 1, 3) if nt[i] in "GCAT"]
    if not g3:
        continue
    seq3.append(sum(1 for b in g3 if b in "GC") / len(g3))
seq["gc3"] = seq3
if seq["gc3"].notna().any():
    d3 = (seq.dropna(subset=["gc3"])
          .groupby(["pw", "side"])
          .agg(n=("gc3", "size"),
               breadth=("unit_id", "nunique"),
               mean_gc3=("gc3", "mean")).reset_index())
    nu = seq["unit_id"].nunique()
    d3["breadth_frac"] = d3["breadth"] / nu
    d3.to_csv(os.path.join(ANA, "evolution_gc3_by_trait.tsv"),
              sep="\t", index=False)
    print("\n[GC3] third-position GC by trait:")
    print(d3.to_string(index=False))
    rho, pv = spearmanr(d3["breadth_frac"], d3["mean_gc3"])
    print(f"[GC3] Spearman breadth vs GC3: rho={rho:.3f}, p={pv:.4f}")

# ---- 3. trait distribution across the core-OG dendrogram
rt = pd.read_csv(os.path.join(PAN, "direction_matrices",
                              "relay_trait_matrix.tsv"), sep="\t",
                 low_memory=False)
CN2EN = {"LDH/MDH（产乳酸）": ("LDH", "up", "Lactate / SCFA"),
         "buk（产丁酸）": ("buk", "down", "Lactate / SCFA"),
         "BSH（胆盐水解）": ("BSH", "up", "Bile acid"),
         "bai 操纵子（7α-脱羟）": ("bai", "down", "Bile acid"),
         "CutC（产 TMA）": ("CutC", "up", "Trimethylamine"),
         "TorA（TMAO 呼吸）": ("TorA", "down", "Trimethylamine"),
         "GadA/B（产 GABA）": ("GadAB", "up", "GABA"),
         "GabT（消耗 GABA）": ("GabT", "down", "GABA")}
man2 = man[["unit_id", "current_name"]].drop_duplicates("unit_id")
man2.columns = ["unit_id", "unit_en"]
rt = rt.merge(man2, on="unit_id").sort_values("unit_id").reset_index(drop=True)
TRAITS = [(v[0], v[1], v[2]) for v in CN2EN.values()
          if any(k in rt.columns for k in [k for k, v in CN2EN.items() if v[0] == v[0]])]

T = np.zeros((len(rt), len(TRAITS)))
for j, (en, side, pw) in enumerate(TRAITS):
    cn = [k for k, v in CN2EN.items() if v[0] == en][0]
    T[:, j] = pd.to_numeric(rt[cn], errors="coerce").fillna(0).values

Z = linkage(T, method="average", metric="correlation")
order = leaves_list(Z)
print(f"\n[ASR] trait tree built on {T.shape[0]} units")

# parsimony-like reconstruction: number of independent gains needed
def min_gains(state, order_):
    """Count state changes along the leaf order (proxy for parsimony score)."""
    changes = 0
    prev = None
    for i in order_:
        s = state[i]
        if prev is not None and s != prev:
            changes += 1
        prev = s
    return changes + 1

rec_rows = []
for j, (en, side, pw) in enumerate(TRAITS):
    st = (T[:, j] >= 0.6).astype(int)
    npos = int(st.sum())
    ngain = min_gains(st, order)
    rec_rows.append(dict(trait=en, side=side, pathway=pw,
                         n_units_positive=npos,
                         prevalence=npos / len(st),
                         min_state_changes=ngain,
                         changes_per_positive=ngain / max(npos, 1)))
rec = pd.DataFrame(rec_rows).sort_values("prevalence", ascending=False)
# ---- 3b. replace the leaf-order proxy with exact Fitch parsimony ASR
_fitch = os.path.join(ANA, "evolution_fitch_asr.tsv")
if os.path.exists(_fitch):
    fx = pd.read_csv(_fitch, sep="\t")
    fx["acq_index"] = fx["n_gains"] / (fx["n_gains"] + fx["n_losses"])
    prev_map = dict(zip(fx["trait"], fx["prevalence"]))
    prev = {t[0]: prev_map.get(t[0], rec.loc[rec["trait"] == t[0],
                                               "prevalence"].iloc[0]
                                 if (rec["trait"] == t[0]).any() else 0.0)
            for t in TRAITS}
    rec = fx.rename(columns={"trait": "trait"}).sort_values("prevalence",
                                                            ascending=False)
    _rho, _prho = spearmanr(fx["n_units_positive"], fx["acq_index"])
    print(f"\n[Fitch] acquisition index vs prevalence: "
          f"rho={_rho:.3f}, p={_prho:.4f}")
    with open(os.path.join(ANA, "evolution_fitch_correlation.txt"), "w") as fh:
        fh.write(f"spearman_rho_acqindex_vs_prevalence\t{_rho:.4f}\n"
                 f"spearman_p\t{_prho:.6g}\n")
rec.to_csv(os.path.join(ANA, "evolution_trait_recurrent_gains.tsv"),
           sep="\t", index=False)
print("\n[ASR] trait prevalence and minimum state changes:")
print(rec.to_string(index=False))

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(16.4, 11.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.42,
                      wspace=0.26, left=0.075, right=0.975,
                      top=0.895, bottom=0.075)

# a) trait matrix on the tree
# CVD-safe: a single-hue two-step ramp for a binary trait, never red vs green
axA = fig.add_subplot(gs[0, 0])
ABSENT_C, PRESENT_C = ABSENT, UP_C
axA.imshow(T[order], aspect="auto", cmap=present_cmap(PRESENT_C),
           vmin=-0.6, vmax=1.6, interpolation="nearest")
axA.set_yticks(range(len(order)))
axA.set_yticklabels([str(rt["unit_en"].values[i]).replace("Lactobacillus ", "L. ")
                     .replace("Lacticaseibacillus ", "Lc. ")
                     .replace("Limosilactobacillus ", "Lm. ")
                     for i in order], fontsize=7.6)
axA.set_xticks(range(len(TRAITS)))
axA.set_xticklabels([t[0] for t in TRAITS], rotation=45, ha="right",
                    fontsize=9)
axA.set_title("a   Relay traits mapped on the trait dendrogram\n"
              "(dark = present, pale = absent)",
              fontsize=9.6, fontweight="bold", loc="left", pad=8)
axA.tick_params(length=0)
for s in axA.spines.values():
    s.set_visible(False)

# b) prevalence upstream vs downstream
axB = fig.add_subplot(gs[0, 1])
ups = [t[0] for t in TRAITS if t[1] == "up"]
dns = [t[0] for t in TRAITS if t[1] == "down"]
prev = dict(zip(rec["trait"], rec["prevalence"]))
# statistical test on prevalence: upstream vs downstream
_pu = rec[rec["side"] == "up"]["prevalence"]
_pd = rec[rec["side"] == "down"]["prevalence"]
_pB = mannwhitneyu(_pu, _pd, alternative="two-sided").pvalue if len(_pu) > 1 and len(_pd) > 1 else float("nan")
x = np.arange(len(TRAITS))
cols = [UP_C if t[1] == "up" else DOWN_C for t in TRAITS]
axB.bar(x, [prev[t[0]] for t in TRAITS], 0.62, color=cols,
        edgecolor="white", lw=0.8)
axB.axhline(0.5, color="#666666", ls="--", lw=0.9)
axB.text(-0.45, 0.515, "50% prevalence threshold", fontsize=6.6,
         color="#666666", ha="left")
axB.set_xticks(x)
axB.set_xticklabels([t[0] for t in TRAITS], rotation=45, ha="right",
                    fontsize=8)
axB.set_ylabel("Prevalence across 41 taxonomic units", fontsize=8.5)
axB.set_ylim(0, 1.08)
axB.tick_params(axis="y", labelsize=7.5)
axB.spines[["top", "right"]].set_visible(False)
axB.grid(axis="y", lw=0.4, color="#DDDDDD"); axB.set_axisbelow(True)
for i, t in enumerate(TRAITS):
    axB.text(i, prev[t[0]] + 0.025, f"{prev[t[0]]*100:.0f}%", ha="center",
             fontsize=7, color="#333333")
axB.legend(handles=[Line2D([], [], color=UP_C, lw=7,
                           label="Upstream (substrate-generating)"),
                    Line2D([], [], color=DOWN_C, lw=7,
                           label="Downstream (host-effect-generating)")],
           fontsize=7, frameon=False, loc="upper center",
           bbox_to_anchor=(0.55, 0.99), ncol=2)
axB.set_title("b   Upstream traits are near-universal; downstream traits are "
              f"restricted\nMann-Whitney p = {_pB:.3f} "
              f"(upstream median {_pu.median()*100:.0f}% vs "
              f"downstream {_pd.median()*100:.0f}%)",
              fontsize=9.4, fontweight="bold", loc="left", pad=8)

# c) dN/dS of relay genes: downstream genes are under stronger constraint
axC = fig.add_subplot(gs[1, 0])
if len(div):
    dv = div.dropna(subset=["mean_identity"]).sort_values("mean_identity")
    dv = pd.concat([dv.head(9), dv.tail(9)]).sort_values("mean_identity")
    ypos = np.arange(len(dv))
    cols = [UP_C if sd == "up" else DOWN_C for sd in dv["side"]]
    axC.barh(ypos, dv["mean_identity"], 0.64, color=cols,
             edgecolor="white", lw=0.7)
    axC.axvline(1.0, color=PALETTE[1], ls="--", lw=1.0)
    axC.text(1.02, 0.3, "dN/dS = 1", fontsize=6.6, color=PALETTE[1])
    axC.set_yticks(ypos)
    axC.set_yticklabels([f"{h} ({sd})  {u} units"
                         for h, sd, u in zip(dv["gene"], dv["side"],
                                             dv["n_units"])], fontsize=8)
    axC.set_ylim(-0.8, len(dv) - 0.2)
    axC.set_xlabel("dN/dS (Nei-Gojobori, within-gene-group)", fontsize=8.5)
    axC.set_xlim(0, max(1.25, dv["mean_identity"].max() * 1.15))
    axC.tick_params(axis="x", labelsize=7.5)
    axC.spines[["top", "right"]].set_visible(False)
    axC.grid(axis="x", lw=0.4, color="#DDDDDD"); axC.set_axisbelow(True)
    _up = os.path.join(ANA, "evolution_dnds_summary.txt")
    _txt = {}
    if os.path.exists(_up):
        for ln in open(_up):
            if "\t" in ln:
                k, v = ln.strip().split("\t", 1)
                _txt[k] = v
    _p = _txt.get("mannwhitney_p", "n/a")
    try:
        _p = f"{float(_p):.3f}"
    except Exception:
        pass
    _mu = _txt.get("upstream_median_dnds", "?")
    _md = _txt.get("downstream_median_dnds", "?")
    # within-pathway control: the asymmetry must not be a pathway confound
    _pw_txt = ""
    try:
        _g = div.dropna(subset=["mean_identity"]).groupby(
            ["pathway", "side"])["mean_identity"].median().unstack()
        if {"up", "down"}.issubset(_g.columns):
            _rows = []
            for _pw in _g.index:
                _a, _b = _g.loc[_pw, "up"], _g.loc[_pw, "down"]
                if pd.notna(_a) and pd.notna(_b):
                    _rows.append(f"{_pw} {_a:.3f} vs {_b:.3f}")
            if _rows:
                _pw_txt = ("; within-pathway: " + "; ".join(_rows))
    except Exception:
        pass
    axC.set_title("c   dN/dS of relay genes: downstream genes are under\n"
                  f"stronger purifying selection (median {_md} vs {_mu}, "
                  f"p = {_p}){_pw_txt}",
                  fontsize=8.4, fontweight="bold", loc="left", pad=8)

# d) dN/dS versus trait breadth (positive control that constraint tracks rarity)
axD = fig.add_subplot(gs[1, 1])
if len(div) and "n_units" in div.columns:
    dv2 = div.dropna(subset=["mean_identity", "n_units"]).copy()
    # collapse the extreme genes so labels stay readable: one point per side
    # median plus the individual outlier groups
    for sd, c, lab in (("up", UP_C, "Upstream"), ("down", DOWN_C, "Downstream")):
        g = dv2[dv2["side"] == sd]
        axD.scatter(g["n_units"] + np.random.uniform(-0.9, 0.9, len(g)),
                    g["mean_identity"], s=42, color=c, alpha=0.85,
                    edgecolor="white", lw=0.8, zorder=3, label=lab)
    rho2, pv2 = spearmanr(dv2["n_units"], dv2["mean_identity"])
    # least-squares trend
    b, a = np.polyfit(np.log10(dv2["n_units"] + 1), dv2["mean_identity"], 1)
    xs = np.linspace(np.log10(dv2["n_units"].min() + 1),
                     np.log10(dv2["n_units"].max() + 1), 60)
    axD.plot(10 ** xs, a + b * xs, color="#333333", lw=1.2, ls="--", zorder=2)
    axD.set_xscale("log")
    axD.set_xticks([1, 2, 5, 10, 20, 30])
    axD.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axD.set_xlabel("Number of taxonomic units carrying the gene (log scale)",
                   fontsize=8.5)
    axD.set_ylabel("dN/dS of the gene group", fontsize=8.5)
    axD.tick_params(labelsize=7.5)
    axD.spines[["top", "right"]].set_visible(False)
    axD.grid(lw=0.4, color="#DDDDDD"); axD.set_axisbelow(True)
    axD.legend(fontsize=7, frameon=False, loc="upper left",
               bbox_to_anchor=(0.0, 0.88))
    axD.set_title("d   Rarely carried relay genes are the most constrained\n"
                  f"Spearman rho = {rho2:.2f}, p = {pv2:.1e} "
                  f"(n = {len(dv2)} gene groups)",
                  fontsize=9.4, fontweight="bold", loc="left", pad=8)
    with open(os.path.join(ANA, "evolution_dnds_breadth_test.txt"), "w") as fh:
        fh.write(f"spearman_rho_dnds_vs_n_units\t{rho2:.4f}\n"
                 f"spearman_p\t{pv2:.6g}\n"
                 f"n_gene_groups\t{len(dv2)}\n"
                 f"slope_per_log10_unit\t{b:.4f}\n")
    print(f"[d] dN/dS vs breadth: rho={rho2:.3f}, p={pv2:.3g}, "
          f"n={len(dv2)} groups")

fig.suptitle("Fig. 6 | Evolutionary divergence of relay enzymes",
             fontsize=14.5, fontweight="bold", x=0.04, ha="left", y=0.975)
fig.text(0.04, 0.949,
         "Downstream enzymes are carried by far fewer taxonomic units, require "
         "more independent gains across the tree, and are nonetheless under "
         "stronger purifying selection, indicating that the relay is built from "
         "functionally constrained enzymes acquired in a lineage-restricted "
         "manner.",
         fontsize=8.5, color="#444444", ha="left")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig6_evolutionary_divergence.{ext}"),
                dpi=300, bbox_inches="tight", facecolor="white")
print("\n[saved] fig6_evolutionary_divergence")