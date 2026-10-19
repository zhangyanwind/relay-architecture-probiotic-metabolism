# -*- coding: utf-8 -*-
"""Direction x unit presence matrices (HMM-hit based, clustering-independent).

Outputs (D:/Probiotics_Pangenome/pangenome/direction_matrices/):
  family_gene_labels.tsv   family -> most common gene symbol + hmm name + direction
  unit_family_counts.tsv   unit x family long table (genomes hit / genomes total / genes)
  unit_family_matrix.tsv   unit x family binary presence matrix (fraction of genomes)
  direction_completeness.tsv  per unit x direction: families present, weighted score
"""
import csv
import os
from collections import Counter, defaultdict

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
ANN = os.path.join(ROOT, "annotations")
PG = os.path.join(ROOT, "pangenome")
OUT = os.path.join(PG, "direction_matrices")
os.makedirs(OUT, exist_ok=True)

HITS = os.path.join(ANN, "hmm_hits_bacteria.tsv")
FSTAT = os.path.join(ANN, "family_stats.tsv")
MAN = os.path.join(ROOT, "metadata", "strain_manifest.tsv")

# ---- manifest: genome -> unit
unit_of, unit_name, unit_genomes = {}, {}, defaultdict(set)
with open(MAN, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r.get("group") != "bacteria":
            continue
        unit_of[r["accession"]] = r["unit_id"]
        unit_name[r["unit_id"]] = r["cn_name"]
        unit_genomes[r["unit_id"]].add(r["accession"])

# ---- family metadata
fam_dir, fam_name = {}, {}
with open(FSTAT, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        fam_dir[r["family"]] = r["direction"]
        fam_name[r["family"]] = r["hmm_name"]

# ---- hits: family -> gene labels, unit x family counts
gene_votes = defaultdict(Counter)
unit_fam_genomes = defaultdict(set)          # (unit, family) -> genomes
unit_fam_genes = Counter()                   # (unit, family) -> genes
fam_units_hit = defaultdict(set)             # family -> units
with open(HITS, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        fam = "%s:%s" % (r["source"], r["hmm_acc"])
        g = r["gene"].strip()
        if g and g.upper() != "NA":
            gene_votes[fam][g] += 1
        uid = unit_of.get(r["genome"])
        if uid is None:
            continue
        unit_fam_genomes[(uid, fam)].add(r["genome"])
        unit_fam_genes[(uid, fam)] += 1
        fam_units_hit[fam].add(uid)

# ---- 1. family gene labels
with open(os.path.join(OUT, "family_gene_labels.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["family", "direction", "hmm_name", "gene_label", "n_units_hit",
                "n_genes_total", "gene_votes"])
    for fam in sorted(fam_dir, key=lambda x: (fam_dir[x], x)):
        tot = sum(unit_fam_genes[(u, fam)] for u in unit_genomes)
        top = gene_votes.get(fam, Counter()).most_common(4)
        label = top[0][0] if top else "-"
        w.writerow([fam, fam_dir[fam], fam_name[fam], label,
                    len(fam_units_hit.get(fam, ())), tot,
                    ",".join("%s:%d" % (g, n) for g, n in top)])

# ---- 2. long table
units = sorted(unit_genomes, key=lambda x: int(x))
families = sorted(fam_dir, key=lambda x: (fam_dir[x], -len(fam_units_hit.get(x, ())), x))
with open(os.path.join(OUT, "unit_family_counts.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["unit_id", "cn_name", "direction", "family", "gene_label", "hmm_name",
                "n_genomes_total", "n_genomes_hit", "pct_genomes", "n_genes"])
    for uid in units:
        for fam in families:
            nh = len(unit_fam_genomes.get((uid, fam), ()))
            ng = len(unit_genomes[uid])
            if nh == 0:
                continue
            top = gene_votes.get(fam, Counter()).most_common(1)
            w.writerow([uid, unit_name[uid], fam_dir[fam], fam,
                        top[0][0] if top else "-", fam_name[fam],
                        ng, nh, round(nh / ng, 3), unit_fam_genes[(uid, fam)]])

# ---- 3. binary matrix (unit x family, values = pct of unit genomes hit)
with open(os.path.join(OUT, "unit_family_matrix.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["unit_id", "cn_name", "direction"] + [f.split(":", 1)[1] for f in families])
    for uid in units:
        row = [uid, unit_name[uid]]
        # direction of a unit: mixed by nature; leave blank, family carries it
        row.append("")
        for fam in families:
            nh = len(unit_fam_genomes.get((uid, fam), ()))
            ng = len(unit_genomes[uid])
            row.append(round(nh / ng, 2) if nh else 0)
        w.writerow(row)

# ---- 4. completeness per unit x direction
n_units = len(units)
fam_n_units = {fam: len(fam_units_hit.get(fam, ())) for fam in families}

with open(os.path.join(OUT, "direction_completeness.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["unit_id", "cn_name", "direction", "n_families_total",
                "n_families_present", "n_families_common", "pct_present",
                "weighted_score", "genes_total"])
    for uid in units:
        for d in ("dir1", "dir2", "dir3"):
            fams = [x for x in families if fam_dir[x] == d]
            present = [x for x in fams if (uid, x) in unit_fam_genomes]
            common = [x for x in fams if fam_n_units[x] >= 0.6 * n_units]
            common_present = [x for x in common if (uid, x) in unit_fam_genomes]
            # weighted: gene abundance normalized by family ubiquity (rare families weigh more)
            wsum = sum((n_units - fam_n_units[x] + 1) / n_units for x in present)
            wtot = sum((n_units - fam_n_units[x] + 1) / n_units for x in fams)
            w.writerow([uid, unit_name[uid], d, len(fams), len(present),
                        "%d/%d" % (len(common_present), len(common)),
                        round(len(present) / len(fams), 3) if fams else 0,
                        round(wsum / wtot, 4) if wtot else 0,
                        sum(unit_fam_genes[(uid, x)] for x in fams)])

print("outputs ->", OUT)
print("units:", n_units, "families:", len(families))
for fn in ("family_gene_labels.tsv", "unit_family_counts.tsv",
           "unit_family_matrix.tsv", "direction_completeness.tsv"):
    p = os.path.join(OUT, fn)
    n = sum(1 for _ in open(p, encoding="utf-8")) - 1
    print("  %-32s %d rows" % (fn, n))
