#!/usr/bin/env python3
"""Statistical analysis of relay pathway enzymes in IBDMDB.

Match HUMAnN2 samples with metadata, perform CD vs UC vs nonIBD comparisons.
"""
import csv
import os
import sys
import numpy as np
from scipy import stats

VALIDATION_DIR = r"D:\Probiotics_Pangenome_v5\validation"
ANALYSIS_DIR = os.path.join(VALIDATION_DIR, "analysis")

# Load enzyme abundance data
enzyme_file = os.path.join(ANALYSIS_DIR, "relay_enzyme_abundance_ibdmdb.tsv")
with open(enzyme_file) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    enzyme_names = header[1:]
    samples = []
    enzyme_data = {}
    for row in reader:
        sample_id = row[0]
        samples.append(sample_id)
        enzyme_data[sample_id] = [float(x) for x in row[1:]]

# Load metadata
metadata_file = os.path.join(VALIDATION_DIR, "hmp2_metadata.csv")
sample_diagnosis = {}
sample_data_type = {}
with open(metadata_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        ext_id = row.get("External ID", "")
        diagnosis = row.get("diagnosis", "")
        data_type = row.get("data_type", "")
        if ext_id:
            sample_diagnosis[ext_id] = diagnosis
            sample_data_type[ext_id] = data_type

# Match samples
print("=== Sample Matching ===")
matched = {}
for s in samples:
    # Try direct match
    if s in sample_diagnosis:
        matched[s] = sample_diagnosis[s]
    else:
        # Try with _P suffix
        s_p = s + "_P"
        if s_p in sample_diagnosis:
            matched[s] = sample_diagnosis[s_p]
        else:
            matched[s] = "unknown"

for s in samples:
    print(f"  {s}: {matched[s]} (data_type: {sample_data_type.get(s, sample_data_type.get(s+'_P', 'unknown'))})")

# Group samples
groups = {"CD": [], "UC": [], "nonIBD": []}
for s in samples:
    diag = matched[s]
    if diag in groups:
        groups[diag].append(s)
    elif diag == "nonIBD":
        groups["nonIBD"].append(s)

print(f"\n=== Group Sizes ===")
for g, ss in groups.items():
    print(f"  {g}: {len(ss)} samples")

# Statistical analysis: CD vs UC vs nonIBD
print(f"\n=== Statistical Analysis (Mann-Whitney U test) ===")
print(f"{'Enzyme':<50} {'CD_mean':>12} {'UC_mean':>12} {'nonIBD_mean':>12} {'CD_vs_UC_p':>12} {'CD_vs_nonIBD_p':>14}")

results = []
for i, enz in enumerate(enzyme_names):
    cd_vals = [enzyme_data[s][i] for s in groups["CD"]]
    uc_vals = [enzyme_data[s][i] for s in groups["UC"]]
    nonibd_vals = [enzyme_data[s][i] for s in groups["nonIBD"]]

    cd_mean = np.mean(cd_vals) if cd_vals else 0
    uc_mean = np.mean(uc_vals) if uc_vals else 0
    nonibd_mean = np.mean(nonibd_vals) if nonibd_vals else 0

    # CD vs UC
    if len(cd_vals) >= 3 and len(uc_vals) >= 3:
        try:
            _, p_cd_uc = stats.mannwhitneyu(cd_vals, uc_vals, alternative="two-sided")
        except ValueError:
            p_cd_uc = np.nan
    else:
        p_cd_uc = np.nan

    # CD vs nonIBD
    if len(cd_vals) >= 3 and len(nonibd_vals) >= 3:
        try:
            _, p_cd_nonibd = stats.mannwhitneyu(cd_vals, nonibd_vals, alternative="two-sided")
        except ValueError:
            p_cd_nonibd = np.nan
    else:
        p_cd_nonibd = np.nan

    results.append({
        "enzyme": enz,
        "cd_mean": cd_mean,
        "uc_mean": uc_mean,
        "nonibd_mean": nonibd_mean,
        "p_cd_uc": p_cd_uc,
        "p_cd_nonibd": p_cd_nonibd,
    })

    print(f"{enz:<50} {cd_mean:>12.6f} {uc_mean:>12.6f} {nonibd_mean:>12.6f} {p_cd_uc:>12.4f} {p_cd_nonibd:>14.4f}")

# Save results
output_file = os.path.join(ANALYSIS_DIR, "relay_enzyme_disease_stats_ibdmdb.tsv")
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["enzyme", "CD_mean", "UC_mean", "nonIBD_mean", "p_CD_vs_UC", "p_CD_vs_nonIBD"])
    for r in results:
        writer.writerow([r["enzyme"], r["cd_mean"], r["uc_mean"], r["nonibd_mean"], r["p_cd_uc"], r["p_cd_nonibd"]])
print(f"\nResults saved to: {output_file}")

# Pathway-level analysis
print(f"\n=== Pathway-Level Analysis ===")
pathway_file = os.path.join(ANALYSIS_DIR, "relay_pathway_abundance_ibdmdb.tsv")
with open(pathway_file) as f:
    reader = csv.reader(f, delimiter="\t")
    pw_header = next(reader)
    pw_names = pw_header[1:]
    pw_data = {}
    for row in reader:
        pw_data[row[0]] = [float(x) for x in row[1:]]

print(f"{'Pathway':<20} {'CD_mean':>12} {'UC_mean':>12} {'nonIBD_mean':>12} {'p_CD_vs_UC':>12} {'p_CD_vs_nonIBD':>14}")
for i, pw in enumerate(pw_names):
    cd_vals = [pw_data[s][i] for s in groups["CD"]]
    uc_vals = [pw_data[s][i] for s in groups["UC"]]
    nonibd_vals = [pw_data[s][i] for s in groups["nonIBD"]]

    cd_mean = np.mean(cd_vals) if cd_vals else 0
    uc_mean = np.mean(uc_vals) if uc_vals else 0
    nonibd_mean = np.mean(nonibd_vals) if nonibd_vals else 0

    if len(cd_vals) >= 3 and len(uc_vals) >= 3:
        try:
            _, p_cd_uc = stats.mannwhitneyu(cd_vals, uc_vals, alternative="two-sided")
        except ValueError:
            p_cd_uc = np.nan
    else:
        p_cd_uc = np.nan

    if len(cd_vals) >= 3 and len(nonibd_vals) >= 3:
        try:
            _, p_cd_nonibd = stats.mannwhitneyu(cd_vals, nonibd_vals, alternative="two-sided")
        except ValueError:
            p_cd_nonibd = np.nan
    else:
        p_cd_nonibd = np.nan

    print(f"{pw:<20} {cd_mean:>12.6f} {uc_mean:>12.6f} {nonibd_mean:>12.6f} {p_cd_uc:>12.4f} {p_cd_nonibd:>14.4f}")

print("\nDone!")
