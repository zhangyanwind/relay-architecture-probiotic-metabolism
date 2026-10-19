# Code repository — relay architecture of probiotic metabolism

Analysis code for the manuscript *Incomplete by architecture: the relay
architecture of probiotic metabolism*.

## Environment

Python 3.13.12 with:

```
numpy==2.5.4
pandas==3.0.6
scipy==1.18.1
matplotlib==3.11.2
```

Install with `pip install -r requirements.txt`.

## Directory contract

Every script resolves its inputs from a single base directory that holds:

```
<BASE>/metadata/strain_manifest.tsv          351-genome manifest
<BASE>/genomes/<accession>/                  RefSeq assemblies (download_genomes.py)
<BASE>/annotations/                          Prodigal/Barrnap annotations
<BASE>/pangenome/                            clustering outputs
<BASE>/pangenome/direction_matrices/         trait matrices, guild matrices
<BASE>/qc/genome_stats_by_unit.tsv           per-unit assembly statistics
<BASE>/validation/analysis/                  result tables written by the scripts
<BASE>/figures/main_text/                    figure output
```

Set the base path at the top of each script (or via `PROB_ROOT` where supported).

## Pipeline order

**1. Input assembly**

| Step | Script | Output |
|---|---|---|
| Download 351 RefSeq complete genomes | `download_genomes.py` | `genomes/`, download log |
| Select the panel and define units (95% ANI + manual curation) | `select_panel_v5.py` | `metadata/strain_manifest.tsv` |
| Annotate proteins (Prodigal) and RNA (Barrnap) | external; see Methods | `annotations/` |
| Two-level protein clustering (90%/80% family, 50%/50% orthogroup) | `cluster_proteins.py` | `pangenome/family_table.tsv`, `orthogroup_table.tsv` |
| Core and openness statistics | `orthogroup_pass.py`, `core_orthogroups.py` | `pangenome/*summary.txt`, openness curves |

**2. Direction annotation**

| Step | Script | Output |
|---|---|---|
| Build HMM library from Pfam/TIGRFAMs | `build_reference_hmms.py` | HMM profiles |
| Extract the direction-relevant profile subset | `extract_hmm_subsets.py` | profile subset |
| Run hmmsearch and record hits | `annotate_pathways.py` | `annotations/hmm_hits*.tsv` |
| Curate 96 profiles into marker/family/broad classes | `curate_families.py` | `family_curation.tsv` |
| Build the functional guild matrix | `functional_guilds.py` | `unit_functional_guilds.tsv` |
| Build the 8-trait relay matrix (41 units) | `relay_figures.py` | `relay_trait_matrix.tsv` |
| Per-unit direction counts and completeness | `direction_matrices.py`, `orthogroup_pass.py` | `unit_direction_genes_orthogroup.tsv`, `direction_completeness.tsv` |

**3. Analysis**

| Step | Script | Figure / table |
|---|---|---|
| Genome-size control | `analysis/r1_1_genome_size_control.py` | Fig. 3e,f; Supplementary Table 10 |
| Pangenome tree and core heat map | `fig1_fig2_core_tree_heatmap.py` | Fig. 1 |
| Pathway module clustering | `fig2_pathway_cluster_heatmap.py` | Fig. 2 |
| Upstream/downstream asymmetry | `fig3_asymmetry.py` | Fig. 3a–d |
| Partner network from metagenomes | `relay_network_analysis.py`, `fig4_relay_network.py` | Fig. 4 |
| Clinical step-resolved disruption | `analyze_metabolome_disease_ibdmdb.py`, `fig5_relay_disruption.py` | Fig. 5 |
| Path-position versus pathway identity | `fig8_pathway_metabolome_collapse.py` | Fig. 8 |
| Host relay capacity | `fig9_host_relay_capacity.py` | Fig. 9 |
| Fitch ancestral states | `fitch_asr_analysis.py` | Fig. 6a,b |
| dN/dS under purifying selection | `evolution_dnds_analysis.py`, `fig6_evolutionary_divergence.py` | Fig. 6c,d |
| Relay completion and orphan carriers | `fig7_relay_completion_analysis.py`, `fig7_first_runner.py` | Fig. 7 |
| Supplementary tables 7–11 | `analysis/make_supp_tables_7_11.py` | Supplementary Tables |
| Supplementary figures | `generate_supp_figures.py` | Supplementary Figs. 1–8 |

## Statistical conventions

- Trait positivity: a unit counts as carrying a trait when at least 60% of its
  genomes carry the gene. This threshold is applied in *every* frequency,
  quadrant and completion analysis, including the Fitch reconstruction.
- Two-level clustering: within-unit analyses use the 90% family layer; all
  cross-unit comparisons use the 50% orthogroup layer.
- Group comparisons: two-sided Mann–Whitney U; paired per-metagenome
  comparisons use Wilcoxon signed-rank; Benjamini–Hochberg within each
  analysis family.
- Colour palette: fixed 13-colour scheme in `paper_palette.py`, with semantic
  colour threading (upstream orange, downstream red, partner green, tryptophan
  purple, bile acid blue, lactate/SCFA yellow, GABA light green, TMA/TMAO teal).

## Licence

MIT (see `LICENSE`).
