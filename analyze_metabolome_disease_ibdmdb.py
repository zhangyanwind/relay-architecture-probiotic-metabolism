#!/usr/bin/env python3
"""Metabolome-disease analysis for IBDMDB (v2 - uses observation metadata).

Load the biom metabolomics file, extract key relay metabolites using
observation metadata (Metabolite name field), compare across CD/UC/nonIBD,
and correlate with enzyme abundance.
"""
import csv
import os
import sys
import numpy as np
from scipy import stats

VALIDATION_DIR = r"D:\Probiotics_Pangenome_v5\validation"
ANALYSIS_DIR = os.path.join(VALIDATION_DIR, "analysis")

# Key metabolites for relay pathways
# Format: (search_term, pathway, role)
# role: "intermediate", "end_product", or "substrate"
KEY_METABOLITES = [
    # Lactate -> SCFA
    ("lactate", "lactate_scfa", "intermediate"),
    ("butyrate", "lactate_scfa", "end_product"),
    ("propionate", "lactate_scfa", "end_product"),
    ("acetate", "lactate_scfa", "end_product"),
    ("caproate", "lactate_scfa", "end_product"),
    ("valerate", "lactate_scfa", "end_product"),
    ("isobutyrate", "lactate_scfa", "end_product"),
    ("isovalerate", "lactate_scfa", "end_product"),
    # Bile acids
    ("cholate", "bile_acid", "intermediate"),
    ("chenodeoxycholate", "bile_acid", "intermediate"),
    ("deoxycholate", "bile_acid", "end_product"),
    ("lithocholate", "bile_acid", "end_product"),
    ("ursodeoxycholate", "bile_acid", "end_product"),
    ("glycocholate", "bile_acid", "intermediate"),
    ("taurocholate", "bile_acid", "intermediate"),
    ("glycodeoxycholate", "bile_acid", "end_product"),
    ("taurodeoxycholate", "bile_acid", "end_product"),
    ("glycochenodeoxycholate", "bile_acid", "intermediate"),
    ("taurochenodeoxycholate", "bile_acid", "intermediate"),
    # TMA/TMAO
    ("trimethylamine-N-oxide", "tma_tmao", "end_product"),
    ("choline", "tma_tmao", "substrate"),
    ("betaine", "tma_tmao", "substrate"),
    ("carnitine", "tma_tmao", "substrate"),
    ("butyrobetaine", "tma_tmao", "substrate"),
    ("gamma-butyrobetaine", "tma_tmao", "substrate"),
    # GABA/glutamate
    ("aminoisobutyric acid/GABA", "gaba_glutamate", "intermediate"),
    ("glutamate", "gaba_glutamate", "substrate"),
    ("glutamine", "gaba_glutamate", "substrate"),
    ("succinate", "gaba_glutamate", "end_product"),
    ("succinate semialdehyde", "gaba_glutamate", "end_product"),
    # Tryptophan
    ("tryptophan", "tryptophan", "substrate"),
    ("indoleacetate", "tryptophan", "end_product"),
    ("indole-3-propionate", "tryptophan", "end_product"),
    ("indole-3-lactate", "tryptophan", "end_product"),
    ("indole-3-acrylate", "tryptophan", "end_product"),
    ("indole", "tryptophan", "end_product"),
    ("kynurenine", "tryptophan", "end_product"),
    ("kynurenate", "tryptophan", "end_product"),
    ("serotonin", "tryptophan", "end_product"),
    ("5-hydroxyindoleacetate", "tryptophan", "end_product"),
    ("5-hydroxytryptophan", "tryptophan", "end_product"),
]


def load_biom_with_metadata(biom_path):
    """Load biom file and return metabolite info with names."""
    import biom
    table = biom.load_table(biom_path)

    metabolites = list(table.ids(axis="observation"))
    samples = list(table.ids(axis="sample"))
    matrix = table.matrix_data.toarray()

    # Build metabolite name mapping from observation metadata
    met_names = {}
    for met_id in metabolites:
        md = table.metadata(met_id, axis="observation")
        if md and "Metabolite" in md:
            met_names[met_id] = md["Metabolite"]
        else:
            met_names[met_id] = met_id

    return metabolites, samples, matrix, met_names


def load_metadata():
    """Load sample metadata."""
    metadata_file = os.path.join(VALIDATION_DIR, "hmp2_metadata.csv")
    sample_info = {}
    with open(metadata_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ext_id = row.get("External ID", "")
            if ext_id:
                sample_info[ext_id] = {
                    "diagnosis": row.get("diagnosis", ""),
                    "data_type": row.get("data_type", ""),
                }
    return sample_info


def load_enzyme_abundance():
    """Load enzyme abundance data from HUMAnN2 analysis."""
    enzyme_file = os.path.join(ANALYSIS_DIR, "relay_enzyme_abundance_ibdmdb.tsv")
    if not os.path.exists(enzyme_file):
        print(f"Warning: enzyme file not found: {enzyme_file}")
        return None, None, None

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

    return enzyme_names, samples, enzyme_data


def main():
    biom_file = os.path.join(VALIDATION_DIR, "HMP2_metabolomics_w_metadata.biom")

    if not os.path.exists(biom_file):
        print(f"ERROR: {biom_file} not found.")
        sys.exit(1)

    # Load data
    print("Loading biom metabolomics data with metadata...")
    metabolites, mbx_samples, matrix, met_names = load_biom_with_metadata(biom_file)
    print(f"  Metabolites: {len(metabolites)}")
    print(f"  Samples: {len(mbx_samples)}")

    # Load metadata
    print("Loading sample metadata...")
    sample_info = load_metadata()
    print(f"  Metadata entries: {len(sample_info)}")

    # Match samples to diagnosis
    sample_diagnosis = {}
    for s in mbx_samples:
        if s in sample_info:
            sample_diagnosis[s] = sample_info[s]["diagnosis"]
        elif s + "_P" in sample_info:
            sample_diagnosis[s] = sample_info[s + "_P"]["diagnosis"]
        else:
            sample_diagnosis[s] = "unknown"

    # Group samples
    groups = {"CD": [], "UC": [], "nonIBD": []}
    for s in mbx_samples:
        diag = sample_diagnosis.get(s, "unknown")
        if diag in groups:
            groups[diag].append(s)

    print(f"\n=== Group Sizes ===")
    for g, ss in groups.items():
        print(f"  {g}: {len(ss)}")

    # Find key metabolites using metadata names
    print("\n=== Key Metabolite Detection (via metadata) ===")
    found_metabolites = {}  # met_id -> (display_name, pathway, role)

    for met_id in metabolites:
        name = met_names.get(met_id, met_id).lower()
        for search_term, pathway, role in KEY_METABOLITES:
            if search_term.lower() == name:  # exact match
                found_metabolites[met_id] = (met_names[met_id], pathway, role)
                break

    print(f"  Found {len(found_metabolites)} key metabolites (exact match):")
    for met_id, (name, pw, role) in sorted(found_metabolites.items(), key=lambda x: x[1][1]):
        print(f"    {name:<40} [{met_id}] -> {pw} ({role})")

    # Statistical analysis: metabolite levels across groups
    if found_metabolites:
        print(f"\n=== Metabolite-Disease Analysis (Mann-Whitney U) ===")
        print(f"{'Metabolite':<40} {'Pathway':<15} {'Role':<12} {'CD_mean':>10} {'UC_mean':>10} {'nonIBD_mean':>12} {'p_CD_UC':>10} {'p_CD_nonIBD':>12} {'p_UC_nonIBD':>12}")

        results = []
        for met_id, (name, pathway, role) in sorted(found_metabolites.items(), key=lambda x: (x[1][1], x[1][0])):
            met_idx = metabolites.index(met_id)

            cd_vals = [matrix[met_idx, mbx_samples.index(s)] for s in groups["CD"] if s in mbx_samples]
            uc_vals = [matrix[met_idx, mbx_samples.index(s)] for s in groups["UC"] if s in mbx_samples]
            nonibd_vals = [matrix[met_idx, mbx_samples.index(s)] for s in groups["nonIBD"] if s in mbx_samples]

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

            # UC vs nonIBD
            if len(uc_vals) >= 3 and len(nonibd_vals) >= 3:
                try:
                    _, p_uc_nonibd = stats.mannwhitneyu(uc_vals, nonibd_vals, alternative="two-sided")
                except ValueError:
                    p_uc_nonibd = np.nan
            else:
                p_uc_nonibd = np.nan

            results.append({
                "metabolite": name,
                "met_id": met_id,
                "pathway": pathway,
                "role": role,
                "cd_mean": cd_mean,
                "uc_mean": uc_mean,
                "nonibd_mean": nonibd_mean,
                "p_cd_uc": p_cd_uc,
                "p_cd_nonibd": p_cd_nonibd,
                "p_uc_nonibd": p_uc_nonibd,
            })

            sig_cu = "***" if p_cd_uc < 0.001 else "**" if p_cd_uc < 0.01 else "*" if p_cd_uc < 0.05 else ""
            sig_cn = "***" if p_cd_nonibd < 0.001 else "**" if p_cd_nonibd < 0.01 else "*" if p_cd_nonibd < 0.05 else ""
            sig_un = "***" if p_uc_nonibd < 0.001 else "**" if p_uc_nonibd < 0.01 else "*" if p_uc_nonibd < 0.05 else ""
            print(f"{name:<40} {pathway:<15} {role:<12} {cd_mean:>10.2f} {uc_mean:>10.2f} {nonibd_mean:>12.2f} {p_cd_uc:>10.4f}{sig_cu:<3} {p_cd_nonibd:>12.4f}{sig_cn:<3} {p_uc_nonibd:>12.4f}{sig_un}")

        # Save results
        output_file = os.path.join(ANALYSIS_DIR, "relay_metabolite_disease_stats_ibdmdb.tsv")
        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["metabolite", "met_id", "pathway", "role", "CD_mean", "UC_mean", "nonIBD_mean", "p_CD_vs_UC", "p_CD_vs_nonIBD", "p_UC_vs_nonIBD"])
            for r in results:
                writer.writerow([r["metabolite"], r["met_id"], r["pathway"], r["role"], r["cd_mean"], r["uc_mean"], r["nonibd_mean"], r["p_cd_uc"], r["p_cd_nonibd"], r["p_uc_nonibd"]])
        print(f"\nResults saved to: {output_file}")

    # Enzyme-metabolite coupling analysis (L1)
    print("\n=== Enzyme-Metabolite Coupling (L1) ===")
    enzyme_names, enz_samples, enzyme_data = load_enzyme_abundance()

    if enzyme_names and found_metabolites:
        # Find overlapping samples (MGX ∩ MBX)
        overlap = set(enz_samples) & set(mbx_samples)
        print(f"  Overlapping MGX+MBX samples: {len(overlap)}")

        if len(overlap) >= 5:
            overlap_list = sorted(overlap)
            print(f"  Samples: {overlap_list}")

            coupling_results = []

            # Define enzyme-metabolite pairs for coupling analysis
            # (enzyme_EC, metabolite_search, expected_correlation, pathway)
            coupling_pairs = [
                ("1.1.1.27", "lactate", "positive", "lactate_scfa"),      # LDH vs lactate
                ("2.7.2.7", "butyrate", "positive", "lactate_scfa"),       # buk vs butyrate
                ("3.5.1.24", "cholate", "positive", "bile_acid"),          # BSH vs cholate
                ("1.7.2.2", "trimethylamine-N-oxide", "positive", "tma_tmao"),  # TorA vs TMAO
                ("4.1.1.15", "aminoisobutyric acid/GABA", "positive", "gaba_glutamate"),  # GadA/B vs GABA
                ("2.6.1.19", "succinate", "positive", "gaba_glutamate"),   # GabT vs succinate
                ("4.2.1.20", "tryptophan", "positive", "tryptophan"),      # trp synthase vs tryptophan
            ]

            for enz_ec, met_search, expected, pathway in coupling_pairs:
                # Find matching enzyme by EC number
                enz_idx = None
                enz_name = None
                for i, en in enumerate(enzyme_names):
                    if enz_ec in en:
                        enz_idx = i
                        enz_name = en
                        break

                # Find matching metabolite
                met_idx = None
                met_display = None
                for met_id, (name, pw, role) in found_metabolites.items():
                    if met_search.lower() == name.lower():
                        met_idx = metabolites.index(met_id)
                        met_display = name
                        break

                if enz_idx is not None and met_idx is not None:
                    enz_vals = [enzyme_data[s][enz_idx] for s in overlap_list]
                    met_vals = [matrix[met_idx, mbx_samples.index(s)] for s in overlap_list if s in mbx_samples]

                    if len(enz_vals) == len(met_vals) and len(enz_vals) >= 5:
                        try:
                            rho, p = stats.spearmanr(enz_vals, met_vals)
                            coupling_results.append({
                                "enzyme": enz_name,
                                "enzyme_ec": enz_ec,
                                "metabolite": met_display,
                                "pathway": pathway,
                                "n": len(enz_vals),
                                "spearman_rho": rho,
                                "p_value": p,
                                "expected": expected,
                            })
                            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                            print(f"  {enz_name:<30} vs {met_display:<30} rho={rho:>7.3f} p={p:>10.4f} {sig}")
                        except Exception as e:
                            print(f"  Error correlating {enz_name} vs {met_display}: {e}")
                    else:
                        print(f"  Skipped {enz_name} vs {met_display}: insufficient data ({len(enz_vals)} values)")
                else:
                    if enz_idx is None:
                        print(f"  Enzyme {enz_ec} not found in HUMAnN2 data")
                    if met_idx is None:
                        print(f"  Metabolite '{met_search}' not found in biom data")

            if coupling_results:
                output_file = os.path.join(ANALYSIS_DIR, "relay_enzyme_metabolite_coupling_ibdmdb.tsv")
                with open(output_file, "w", newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerow(["enzyme", "enzyme_ec", "metabolite", "pathway", "n_samples", "spearman_rho", "p_value", "expected_direction"])
                    for r in coupling_results:
                        writer.writerow([r["enzyme"], r["enzyme_ec"], r["metabolite"], r["pathway"], r["n"], r["spearman_rho"], r["p_value"], r["expected"]])
                print(f"\n  Coupling results saved to: {output_file}")
            else:
                print("  No coupling pairs could be analyzed")
        else:
            print(f"  Insufficient overlapping samples ({len(overlap)} < 5) for coupling analysis")
    else:
        print("  Enzyme data or metabolites not available for coupling analysis")

    print("\nDone!")


if __name__ == "__main__":
    main()
