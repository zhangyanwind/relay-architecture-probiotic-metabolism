#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate Supplementary Tables 7-11 for the second-round revision.

  7  gene -> profile -> pathway mapping with enzyme commission numbers
  8  per-unit genome assembly statistics (source of the Fig. 3e,f control)
  9  direction curation table: family, evidence class, marker-gene share
 10  genome-size control analysis, per unit
 11  Fig. 8 / Fig. 9 source values in one place
"""
import os
import shutil
import numpy as np
import pandas as pd

ROOT = "D:/Probiotics_Pangenome_v5"
DM = f"{ROOT}/pangenome/direction_matrices"
ANA = f"{ROOT}/validation/analysis"
QC = f"{ROOT}/qc"
OUT = "C:/Users/zhangyan/WorkBuddy/2026-10-17-17-38-14/npj投稿材料/supplementary_tables"
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------ Supp Table 7
lab = pd.read_csv(f"{DM}/family_gene_labels.tsv", sep="\t")
cur = pd.read_csv(f"{DM}/family_curation.tsv", sep="\t")
# attach an enzyme-commission number where the curated label defines one
EC = {
    "Ldh_1_N": "EC 1.1.1.27", "But_2": "EC 2.7.2.7", "BSH": "EC 3.5.1.26",
    "BaiE": "EC 1.6.99.7", "CutC": "EC 1.14.13.164", "TorA": "EC 1.6.99.9",
    "GadA": "EC 4.1.1.15", "GabT": "EC 2.6.1.19", "TnaA": "EC 4.1.99.03",
}
PW = {"dir1": "upstream (substrate-generating)",
      "dir2": "downstream (host-effect-generating)",
      "dir3": "other (outside the five profiled pathways)",
      "other": "other"}
t7 = lab.merge(
    cur[["family", "class", "marker_share", "n_hits", "n_named"]],
    on="family", how="left")
t7["enzyme_commission"] = t7.gene_label.map(EC)
t7["metabolic_direction"] = t7.direction.map(PW)
t7 = t7.rename(columns={
    "family": "hmm_profile", "hmm_name": "model_name",
    "gene_label": "curated_gene_label", "n_units_hit": "n_taxonomic_units",
    "n_genes_total": "n_genes_all_genomes", "gene_votes": "gene_symbol_votes",
    "class": "evidence_class", "marker_share": "marker_gene_share"})
t7 = t7[["hmm_profile", "model_name", "curated_gene_label",
         "enzyme_commission", "metabolic_direction", "evidence_class",
         "marker_gene_share", "n_taxonomic_units", "n_genes_all_genomes",
         "n_hits", "n_named", "gene_symbol_votes"]]
t7 = t7.fillna("")
t7.to_csv(f"{OUT}/supp_table7_gene_profile_ec_mapping.csv", index=False)
print(f"[7]  {len(t7)} profiles")

# ------------------------------------------------------------ Supp Table 8
gs = pd.read_csv(f"{QC}/genome_stats_by_unit.tsv", sep="\t")
man = pd.read_csv(f"{ROOT}/metadata/strain_manifest.tsv", sep="\t")
namecol = "current_name" if "current_name" in man.columns else "cn_name"
ud = pd.read_csv(f"{ROOT}/pangenome/unit_direction_genes_orthogroup.tsv", sep="\t")
tot = ud[["dir1_genes", "dir2_genes", "dir3_genes"]].sum(axis=1)
ud["frac_up"] = ud.dir1_genes / tot
ud["frac_down"] = ud.dir2_genes / tot
ud["ratio_up_down"] = (ud.dir1_genes + 1) / (ud.dir2_genes + 1)
t8 = gs.merge(ud[["unit_id", "frac_up", "frac_down", "ratio_up_down"]],
              on="unit_id", how="left")
latin = man.groupby("unit_id")[namecol].first() if "unit_id" in man.columns else None
if latin is not None:
    t8["unit_name_latin"] = t8.unit_id.map(latin)
keep = ["unit_id", "unit_name_latin", "cn_name", "n_genomes", "median_size_mb",
        "min_size_mb", "max_size_mb", "median_gc", "median_n50_mb",
        "median_contigs", "min_pair_ani", "max_pair_ani", "ani_span_pt",
        "n_type_strain", "n_reference", "flag",
        "dir1_genes", "dir2_genes", "frac_up", "frac_down", "ratio_up_down"]
t8 = t8[[c for c in keep if c in t8.columns]]
t8 = t8.rename(columns={"cn_name": "unit_name_chinese",
                        "median_size_mb": "median_assembly_mb",
                        "min_size_mb": "min_assembly_mb",
                        "max_size_mb": "max_assembly_mb",
                        "median_n50_mb": "median_N50_mb",
                        "frac_up": "fraction_genes_upstream",
                        "frac_down": "fraction_genes_downstream",
                        "ratio_up_down": "upstream_downstream_ratio"})
t8.to_csv(f"{OUT}/supp_table8_unit_genome_stats.csv", index=False)
print(f"[8]  {len(t8)} units")

# ------------------------------------------------------------ Supp Table 9
t9 = cur.rename(columns={
    "family": "hmm_profile", "hmm_name": "model_name",
    "class": "evidence_class", "marker_genes": "expected_marker_genes",
    "n_hits": "n_protein_hits", "n_named": "n_hits_with_gene_symbol",
    "marker_share": "marker_gene_share", "top_genes": "most_frequent_gene_symbols"})
t9["metabolic_direction"] = t9.direction.map(PW)
t9["n_protein_hits"] = t9.n_protein_hits.fillna(0).astype(int)
t9["n_hits_with_gene_symbol"] = t9.n_hits_with_gene_symbol.fillna(0).astype(int)
t9["marker_gene_share"] = t9.marker_gene_share.astype(float).round(4)
t9 = t9[["hmm_profile", "model_name", "metabolic_direction", "evidence_class",
         "expected_marker_genes", "n_protein_hits", "n_hits_with_gene_symbol",
         "marker_gene_share", "most_frequent_gene_symbols"]]
t9.to_csv(f"{OUT}/supp_table9_direction_curation.csv", index=False)
print(f"[9]  {len(t9)} curated profiles; "
      f"{int((t9.marker_gene_share >= 0.8).sum())} with marker share >= 0.8")

# ------------------------------------------------------------ Supp Table 10
t10 = pd.read_csv(f"{ANA}/r1_1_genome_size_control.tsv", sep="\t")
t10 = t10.rename(columns={"cn_name": "unit_name_chinese"})
if latin is not None:
    t10["unit_name_latin"] = t10.unit_id.map(latin)
t10["upstream_downstream_ratio"] = t10.ratio_ud
t10.to_csv(f"{OUT}/supp_table10_genome_size_control.csv", index=False)
print(f"[10] {len(t10)} units, genome size "
      f"{t10.genome_mb.min():.2f}-{t10.genome_mb.max():.2f} Mb")

# ------------------------------------------------------------ Supp Table 11
f8 = pd.read_csv(f"{ANA}/fig8_metabolite_lfc_by_pathway.tsv", sep="\t")
f8p = pd.read_csv(f"{ANA}/fig8_pathway_effect_summary.tsv", sep="\t")
f8s = pd.read_csv(f"{ANA}/fig8_step_resolved_summary.tsv", sep="\t")
f9m = pd.read_csv(f"{ANA}/fig9_host_relay_capacity_per_metagenome.tsv", sep="\t")
f9p = pd.read_csv(f"{ANA}/fig9_pathway_capacity_balance.tsv", sep="\t")
t11a = f8.rename(columns={"pathway": "pathway", "name": "pathway_name",
                          "role": "relay_position"})
t11a.to_csv(f"{OUT}/supp_table11a_fig8_metabolite_lfc.csv", index=False)
t11b = f8p[["pathway", "name", "n", "median", "q1", "q3"]].rename(
    columns={"name": "pathway_name", "median": "median_abs_log2FC",
             "q1": "Q1_abs_log2FC", "q3": "Q3_abs_log2FC"})
t11b.to_csv(f"{OUT}/supp_table11b_fig8_pathway_effects.csv", index=False)
t11c = f8s.rename(columns={"step": "relay_position", "cd_med": "median_log2FC_CD",
                          "uc_med": "median_log2FC_UC",
                          "cd_abs": "median_abs_log2FC_CD"})
t11c.to_csv(f"{OUT}/supp_table11c_fig8_step_resolved.csv", index=False)
t11d = f9m.rename(columns={c: c for c in f9m.columns})
t11d.to_csv(f"{OUT}/supp_table11d_fig9_capacity_per_metagenome.csv", index=False)
t11e = f9p.rename(columns={"pathway": "pathway", "name": "pathway_name",
                          "up_med": "median_upstream_abundance",
                          "dn_med": "median_downstream_abundance",
                          "frac_hosts_up_gt_down":
                              "fraction_metagenomes_upstream_exceeds_downstream"})
t11e.to_csv(f"{OUT}/supp_table11e_fig9_pathway_balance.csv", index=False)
print(f"[11] Fig8 {len(t11a)} metabolites / {len(t11b)} pathways / "
      f"{len(t11c)} positions; Fig9 {len(t11d)} metagenomes / {len(t11e)} pathways")

print("\n[done] supplementary tables written to", OUT)
for f in sorted(os.listdir(OUT)):
    print("   ", f)