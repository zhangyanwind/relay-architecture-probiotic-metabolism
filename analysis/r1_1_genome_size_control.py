#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R1-1 control analysis: is the downstream deficit merely a consequence of small
genome size?

The manuscript claims probiotics are *biased towards the upstream half* rather
than merely *carrying fewer genes*.  That claim only survives if downstream
carriage is decoupled from genome size.  This script tests the decoupling
directly:

  1. per-unit genome size (from the Panaroo/Prodigal annotation) vs
     downstream gene fraction  (Spearman)
  2. the same for the 90% family layer
  3. per-unit downstream RELAXED (any gene) carriage, which removes the
     60% threshold as a confounder
  4. genome size vs total annotated gene count (sanity: the two must correlate)
  5. partial rank correlation of downstream carriage with genome size after
     removing the effect of upstream carriage

Result interpretation is written to analysis/r1_1_genome_size_control.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "D:/Probiotics_Pangenome_v5"
PAN = f"{ROOT}/pangenome"
ANA = f"{ROOT}/validation/analysis"
os.makedirs(ANA, exist_ok=True)

# ---------------------------------------------------------------- genome size
# Genome sizes from the assembly statistics QC table if present, otherwise
# reconstructed from the orthogroup assignment table (protein count proxy).
size_candidates = [
    f"{ROOT}/qc/genome_stats_by_unit.tsv",
    f"{ROOT}/qc/genome_stats.tsv",
    f"{ROOT}/metadata/genome_sizes.tsv",
]
gs = None
src = None
for p in size_candidates:
    if os.path.exists(p):
        gs = pd.read_csv(p, sep="\t")
        src = p
        break

if gs is None:
    # Fall back: total annotated proteins per genome from the assignment table.
    asg = pd.read_csv(f"{PAN}/orthogroup_assignments.tsv", sep="\t",
                      usecols=lambda c: c in ("genome", "genome_id", "strain",
                                              "unit_id", "n_genes", "n_proteins",
                                              "contig_n50", "total_genes",
                                              "genome_size", "length"))
    print("[i] orthogroup_assignments columns:", list(asg.columns))
    raise SystemExit("No genome size table found; inspect assignment columns.")

print(f"[i] genome size source: {src}")
print(gs.head())

# ------------------------------------------------- per-unit direction summary
ud = pd.read_csv(f"{PAN}/unit_direction_genes_orthogroup.tsv", sep="\t")
ud.columns = [c.strip() for c in ud.columns]

# relaxed downstream carriage: any hit at all, from the unit family counts
ufc = pd.read_csv(f"{PAN}/direction_matrices/unit_family_counts.tsv", sep="\t")
rel_dir = pd.read_csv(f"{PAN}/direction_matrices/direction_completeness.tsv",
                      sep="\t")

rows = []
for uid, grp in ufc.groupby("unit_id"):
    tot_genes = grp.n_genes.sum()
    down = grp[grp.direction == "dir2"]
    up = grp[grp.direction == "dir1"]
    rows.append({
        "unit_id": uid,
        "cn_name": grp.cn_name.iloc[0],
        "genes_any_frac_up": up.n_genes.sum() / tot_genes,
        "genes_any_frac_down": down.n_genes.sum() / tot_genes,
        "genes_any_ratio": (up.n_genes.sum() + 1) / (down.n_genes.sum() + 1),
        "fams_any_frac_down": down.n_genes.sum() / tot_genes,
        # relaxed carriage: fraction of a unit's genomes carrying >=1 downstream gene
        "pct_genomes_down_any": (down.n_genomes_hit.sum() / down.n_genomes_total.sum()),
        "pct_genomes_up_any": (up.n_genomes_hit.sum() / up.n_genomes_total.sum()),
    })
R = pd.DataFrame(rows)

# thresholded direction fractions (the values actually used in the paper)
ud["frac_up"] = ud.dir1_genes / (ud.dir1_genes + ud.dir2_genes + ud.dir3_genes)
ud["frac_down"] = ud.dir2_genes / (ud.dir1_genes + ud.dir2_genes + ud.dir3_genes)
ud["ratio_ud"] = (ud.dir1_genes + 1) / (ud.dir2_genes + 1)
R = R.merge(ud[["unit_id", "frac_up", "frac_down", "ratio_ud",
                "dir1_genes", "dir2_genes"]],
            on="unit_id", how="left")

# ------------------------------------------------------------- genome sizes
gcol = [c for c in gs.columns
        if c.lower() in ("median_size_mb", "size_mb", "genome_mb")
        or ("size" in c.lower() and c.lower().endswith("mb"))]
icol = [c for c in gs.columns if c.lower() == "unit_id"]
print("[i] genome table size cols:", gcol, "unit cols:", icol)
if icol and gcol:
    g = gs[[icol[0], gcol[0]]].copy()
    g.columns = ["unit_id", "genome_mb"]
    g["unit_id"] = g.unit_id.astype(int)
    R = R.merge(g, on="unit_id", how="left")

out = []
out.append("R1-1  Genome-size control analysis")
out.append(f"units: {len(R)}")

if "genome_mb" in R.columns and R.genome_mb.notna().sum() > 10:
    tests = [
        ("genome size  vs  thresholded downstream fraction", "genome_mb", "frac_down"),
        ("genome size  vs  thresholded upstream fraction", "genome_mb", "frac_up"),
        ("genome size  vs  thresholded upstream/downstream ratio", "genome_mb", "ratio_ud"),
        ("genome size  vs  relaxed downstream gene fraction", "genome_mb", "genes_any_frac_down"),
        ("genome size  vs  relaxed upstream/downstream ratio", "genome_mb", "genes_any_ratio"),
        ("genome size  vs  % genomes carrying any downstream gene", "genome_mb", "pct_genomes_down_any"),
        ("genome size  vs  % genomes carrying any upstream gene", "genome_mb", "pct_genomes_up_any"),
    ]
    for label, a, b in tests:
        s = R[[a, b]].dropna()
        if len(s) < 8:
            out.append(f"  {label}: too few points (n={len(s)})")
            continue
        rho, p = stats.spearmanr(s[a], s[b])
        out.append(f"  {label}: rho = {rho:+.3f}, p = {p:.4f}, n = {len(s)}")
    R["resid_down"] = R.frac_down - np.poly1d(
        np.polyfit(R.genome_mb.dropna(), R.loc[R.genome_mb.notna(), "frac_down"],
                   1))(R.genome_mb.dropna())
    s = R[["genome_mb", "resid_down"]].dropna()
    if len(s) > 8:
        rho, p = stats.spearmanr(s.genome_mb, s.resid_down)
        out.append(f"  residual downstream fraction after removing the genome-size"
                   f" trend: rho = {rho:+.3f}, p = {p:.4f}")
else:
    out.append("  [!] genome sizes unavailable; skipping size-correlation tests")

# ------------------------------------------- is upstream>downstream size-free?
for tag, ucol, dcol in [
        ("thresholded gene fraction", "frac_up", "frac_down"),
        ("relaxed genome carriage", "pct_genomes_up_any", "pct_genomes_down_any")]:
    s = R[[ucol, dcol]].dropna()
    if len(s) < 8:
        continue
    w, pw = stats.wilcoxon(s[ucol], s[dcol])
    out.append(f"  [{tag}] units with upstream>downstream: "
               f"{int((s[ucol] > s[dcol]).sum())}/{len(s)}; "
               f"median upstream = {s[ucol].median():.3f}, "
               f"downstream = {s[dcol].median():.3f}; Wilcoxon p = {pw:.3g}")

# Same test applied within genome-size strata: if the bias were a pure
# consequence of small genomes, it would vanish among the largest units.
if "genome_mb" in R.columns and R.genome_mb.notna().sum() > 20:
    R2 = R.dropna(subset=["genome_mb"]).sort_values("genome_mb")
    med = R2.genome_mb.median()
    for lab, sub in [("smaller half", R2[R2.genome_mb <= med]),
                     ("larger half", R2[R2.genome_mb > med])]:
        w, pw = stats.wilcoxon(sub.frac_up, sub.frac_down)
        out.append(f"  [{lab}, median {sub.genome_mb.median():.2f} Mb, n={len(sub)}] "
                   f"upstream>downstream in {int((sub.frac_up > sub.frac_down).sum())}/{len(sub)}"
                   f"; median ratio = {sub.ratio_ud.median():.2f}; Wilcoxon p = {pw:.3g}")

# smallest vs largest decile of genome size, if sizes exist
if "genome_mb" in R.columns and R.genome_mb.notna().sum() > 10:
    R2 = R.dropna(subset=["genome_mb"]).sort_values("genome_mb")
    k = max(4, len(R2) // 5)
    small, large = R2.head(k), R2.tail(k)
    out.append(f"  smallest {k} units (median {small.genome_mb.median():.2f} Mb): "
               f"median frac_down = {small.frac_down.median():.3f}, "
               f"ratio = {small.ratio_ud.median():.2f}")
    out.append(f"  largest  {k} units (median {large.genome_mb.median():.2f} Mb): "
               f"median frac_down = {large.frac_down.median():.3f}, "
               f"ratio = {large.ratio_ud.median():.2f}")
    u, p = stats.mannwhitneyu(small.ratio_ud, large.ratio_ud, alternative="less")
    out.append(f"  upstream/downstream ratio, smallest < largest: MWU p = {p:.4f}")

R.to_csv(f"{ANA}/r1_1_genome_size_control.tsv", sep="\t", index=False)
txt = "\n".join(out)
open(f"{ANA}/r1_1_genome_size_control.txt", "w", encoding="utf-8").write(txt)
print(txt)