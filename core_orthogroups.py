# -*- coding: utf-8 -*-
"""What is inside the 40 collection-wide core orthogroups?

For each core orthogroup (present in >=99% of the 351 genomes) report:
  * number of proteins, number of genomes
  * the largest 90% family inside it and its gene symbol / product (from the
    GFF annotation of the representative protein)
  * how many 90% families it contains

Output: pangenome/core_orthogroups.tsv  +  console table
"""
import csv
import os
import re
from collections import defaultdict

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome_v5")
PG = os.path.join(ROOT, "pangenome")
ANN = os.path.join(ROOT, "proteins", "annotation_all.tsv")
OUT = os.path.join(PG, "core_orthogroups.tsv")

# representative protein of each 90% family: "ACC|protein_id" -> n_proteins
fam_rep = {}
fam_n = {}
fam_cat = {}
with open(os.path.join(PG, "family_table.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        fam_rep[int(r["family_id"])] = r["representative_protein"]
        fam_n[int(r["family_id"])] = int(r["n_proteins"])
        fam_cat[int(r["family_id"])] = r["category"]

# position join family <-> orthogroup (same row order, exact)
og_of_fam = {}
with open(os.path.join(PG, "gene_family_assignments.tsv"), encoding="utf-8") as f1, \
     open(os.path.join(PG, "orthogroup_assignments.tsv"), encoding="utf-8") as f2:
    for a, b in zip(csv.DictReader(f1, delimiter="\t"), csv.DictReader(f2, delimiter="\t")):
        og_of_fam[int(a["family_id"])] = int(b["orthogroup_id"])

og_info = {}
with open(os.path.join(PG, "orthogroup_table.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        og_info[int(r["orthogroup_id"])] = r

# gene symbol / product lookup for loci (first non-empty wins)
sym = {}
prod = {}
with open(ANN, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["type"] != "CDS":
            continue
        k = (r["accession"], r["locus_tag"])
        if k in sym:
            continue
        if r["gene"]:
            sym[k] = r["gene"]
        if r["product"]:
            prod[k] = r["product"]

og_fams = defaultdict(list)
for fam, og in og_of_fam.items():
    og_fams[og].append(fam)

rows = []
for og, info in og_info.items():
    if info["category"] != "core":
        continue
    fams = sorted(og_fams.get(og, []), key=lambda x: -fam_n.get(x, 0))
    biggest = fams[0] if fams else None
    gname = ""
    pname = ""
    if biggest is not None:
        rep = fam_rep.get(biggest, "")
        parts = rep.split("|")
        if len(parts) >= 2:
            k = (parts[0], parts[1])
            gname = sym.get(k, "")
            pname = prod.get(k, "")
    rows.append({
        "orthogroup_id": og,
        "n_proteins": int(info["n_proteins"]),
        "n_genomes": int(info["n_genomes"]),
        "genome_fraction": info["genome_fraction"],
        "n_families_90": len(fams),
        "biggest_family": biggest,
        "biggest_family_size": fam_n.get(biggest, "") if biggest is not None else "",
        "gene_symbol": gname,
        "product": (pname or "")[:70],
    })
rows.sort(key=lambda r: -r["n_proteins"])

with open(OUT, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t",
                       lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

print("core orthogroups: %d, proteins %d" % (len(rows), sum(r["n_proteins"] for r in rows)))
print("%-6s %6s %5s %5s %-10s %-28s %s" % ("og", "prot", "gen", "fam#", "gene", "product", "fam_size"))
for r in rows:
    print("%-6d %6d %5d %5d %-10s %-28s %s" % (
        r["orthogroup_id"], r["n_proteins"], r["n_genomes"], r["n_families_90"],
        (r["gene_symbol"] or "-")[:10], (r["product"] or "-")[:28], r["biggest_family_size"]))
print("saved ->", OUT)
