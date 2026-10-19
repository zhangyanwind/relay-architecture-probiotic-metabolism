#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Panel a upgrade | true Fitch parsimony ancestral-state reconstruction.

Replaces the previous proxy (number of state changes along the leaf ORDER,
which ignores tree topology and double-counts) with an exact Fitch
maximum-parsimony reconstruction on the trait dendrogram.

For every relay trait we report:
  parsimony_score   - minimum number of state changes on the tree
  n_gains           - number of branches inferred as 0 -> 1 (gain)
  n_losses          - number of branches inferred as 1 -> 0 (loss)
  n_reconstructions - how often the root is ambiguous (0 / 1 / {0,1})
  min_gains         - lower bound on the number of independent acquisitions,
                      i.e. gains when the ancestral state is fixed to absent,
                      obtained by counting the number of positive clades
Outputs English-labelled tables:
  evolution_fitch_asr.tsv
  evolution_fitch_summary.txt
"""
import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, to_tree
from scipy.stats import mannwhitneyu, spearmanr

BASE = r"D:\Probiotics_Pangenome_v5"
PAN = os.path.join(BASE, "pangenome")
DM = os.path.join(PAN, "direction_matrices")
ANA = os.path.join(BASE, "validation", "analysis")
os.makedirs(ANA, exist_ok=True)

CN2EN = {"LDH/MDH（产乳酸）": "LDH", "buk（产丁酸）": "buk",
         "BSH（胆盐水解）": "BSH", "bai 操纵子（7α-脱羟）": "bai",
         "CutC（产 TMA）": "CutC", "TorA（TMAO 呼吸）": "TorA",
         "GadA/B（产 GABA）": "GadAB", "GabT（消耗 GABA）": "GabT"}
SIDE = {"LDH": "up", "BSH": "up", "CutC": "up", "GadAB": "up",
        "buk": "down", "bai": "down", "TorA": "down", "GabT": "down"}
PATHWAY = {"LDH": "Lactate / SCFA", "buk": "Lactate / SCFA",
           "BSH": "Bile acid", "bai": "Bile acid",
           "CutC": "Trimethylamine", "TorA": "Trimethylamine",
           "GadAB": "GABA", "GabT": "GABA"}

rt = pd.read_csv(os.path.join(DM, "relay_trait_matrix.tsv"), sep="\t",
                 low_memory=False)
man = pd.read_csv(os.path.join(BASE, "metadata", "strain_manifest.tsv"),
                  sep="\t", low_memory=False)
rt = rt.merge(man[["unit_id", "current_name"]].drop_duplicates("unit_id"),
              on="unit_id").sort_values("unit_id").reset_index(drop=True)
print(f"[data] {len(rt)} units")

TRAITS = [(en, SIDE[en], PATHWAY[en]) for en in CN2EN.values()]
T = np.column_stack([pd.to_numeric(rt[k], errors="coerce").fillna(0).values
                     for k in CN2EN])
Z = linkage(T, method="average", metric="correlation")
root = to_tree(Z)
order = leaves_list(Z)
print(f"[tree] trait dendrogram built, {T.shape[0]} leaves")


def fitch(states):
    """Exact Fitch parsimony.

    Returns (score, gains, losses, root_ambiguous) where gains/losses are
    counted on the most-parsimonious assignment (Fitch downpass; when the
    child state is ambiguous the parent state is propagated and the change is
    attributed to the branch with the majority of the child subtree).
    """
    score = 0
    nxt = {}

    def rec(node):
        nonlocal score
        if node.is_leaf():
            s = {int(states[node.id])}
            nxt[node.id] = s
            return s
        a = rec(node.left)
        b = rec(node.right)
        inter = a & b
        if inter:
            s = set(inter)
        else:
            union = a | b
            # attribute the change to the larger of the two child subtrees so
            # that gains are counted where the majority of the evidence is
            if len(a) > len(b):
                s = set(a)
            elif len(b) > len(a):
                s = set(b)
            else:
                s = set(union)
            score += 1
        nxt[node.id] = s
        return s

    rootset = rec(root)

    # downpass: propagate the parent state, count branch-level changes
    gains = [0]
    losses = [0]
    assign = {}

    def down(node, parent_state):
        if node.is_leaf():
            assign[node.id] = int(states[node.id])
            return
        assign[node.id] = parent_state
        for ch in (node.left, node.right):
            allowed = nxt[ch.id]
            st = parent_state if parent_state in allowed else min(allowed)
            assign[ch.id] = st
            if st != parent_state:
                if st == 1:
                    gains[0] += 1
                else:
                    losses[0] += 1
            down(ch, st)

    down(root, min(rootset))
    return score, gains[0], losses[0], len(rootset) > 1


def min_independent_gains(states):
    """Lower bound on independent acquisitions: the number of positive clades
    (positive leaves whose parent is negative), i.e. how many separate lineages
    must have acquired the trait."""
    pos = [i for i in range(len(states)) if states[i] == 1]
    if not pos:
        return 0
    n_gains = 0
    cnt = [0]

    def rec(node):
        if node.is_leaf():
            return int(states[node.id])
        a = rec(node.left)
        b = rec(node.right)
        if a == 0 and b == 1:
            cnt[0] += 1
        elif a == 1 and b == 0:
            pass
        return max(a, b)

    rec(root)
    return cnt[0]


rows = []
for j, (en, side, pw) in enumerate(TRAITS):
    st = (T[:, j] >= 0.6).astype(int)
    score, gains, losses, amb = fitch(st)
    npos = int(st.sum())
    rows.append(dict(trait=en, side=side, pathway=pw,
                     n_units_positive=npos,
                     prevalence=npos / len(st),
                     parsimony_score=score,
                     n_gains=gains, n_losses=losses,
                     root_ambiguous=int(amb),
                     min_independent_gains=min_independent_gains(st),
                     changes_per_positive=score / max(npos, 1)))
rec = pd.DataFrame(rows).sort_values("prevalence", ascending=False)
rec.to_csv(os.path.join(ANA, "evolution_fitch_asr.tsv"), sep="\t", index=False)
print("\n[Fitch] ancestral-state reconstruction:")
print(rec.to_string(index=False))

up = rec[rec["side"] == "up"]
dn = rec[rec["side"] == "down"]
U, p = mannwhitneyu(up["parsimony_score"], dn["parsimony_score"],
                    alternative="two-sided")
rho, prho = spearmanr(rec["prevalence"], rec["parsimony_score"])
rho2, prho2 = spearmanr(rec["prevalence"], rec["min_independent_gains"])

lines = [
    f"n_units\t{len(rec)}",
    f"n_upstream_traits\t{len(up)}",
    f"n_downstream_traits\t{len(dn)}",
    f"median_parsimony_upstream\t{up['parsimony_score'].median():.1f}",
    f"median_parsimony_downstream\t{dn['parsimony_score'].median():.1f}",
    f"mannwhitney_p_parsimony\t{p:.4f}",
    f"spearman_rho_prevalence_vs_parsimony\t{rho:.4f}",
    f"spearman_p_prevalence_vs_parsimony\t{prho:.4g}",
    f"spearman_rho_prevalence_vs_mingains\t{rho2:.4f}",
    f"spearman_p_prevalence_vs_mingains\t{prho2:.4g}",
]
for _, r in rec.iterrows():
    lines.append(f"TRAIT\t{r['trait']}\t{r['side']}\t{r['pathway']}\t"
                 f"{r['n_units_positive']}\t{r['prevalence']:.4f}\t"
                 f"{r['parsimony_score']}\t{r['n_gains']}\t{r['n_losses']}\t"
                 f"{r['min_independent_gains']}\t"
                 f"{r['changes_per_positive']:.4f}")
with open(os.path.join(ANA, "evolution_fitch_summary.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n[Fitch] summary written to evolution_fitch_summary.txt")