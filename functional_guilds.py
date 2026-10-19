# -*- coding: utf-8 -*-
"""Functional trait matrix per species unit (paper-ready).

Groups the curated marker families into pathway-level traits, then scores each
unit by the fraction of its genomes carrying the trait's marker families.

Outputs:
  unit_functional_guilds.tsv   unit x trait (fraction of genomes, 0-1)
  unit_guild_summary.txt       per-trait carrier list + guild composition
"""
import csv
import os

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
DM = os.path.join(ROOT, "pangenome", "direction_matrices")

# trait -> (direction, [family short names], rule)
#   rule: "any" (>=1 marker) or "majority" (>=half of the listed markers)
TRAITS = [
    # ---------------- dir1 ----------------
    ("乳酸脱氢酶 LDH/MDH 家族", "dir1", ["PF00056", "TIGR01763", "TIGR01771"], "any"),
    ("D-乳酸脱氢酶", "dir1", ["PF09330"], "any"),
    ("乳酸消旋酶 LarA", "dir1", ["PF09861", "PF21113"], "any"),
    ("乳酸转运蛋白", "dir1", ["PF02652", "TIGR00795", "LCTP"], "any"),
    ("乙酸生成 pta-ackA", "dir1", ["TIGR00651", "TIGR02706", "TIGR00016"], "majority"),
    ("乙醇生成 adhE", "dir1", ["ADHE", "TIGR02518", "TIGR01722"], "any"),
    ("甲酸生成 pflB", "dir1", ["TIGR01255", "PFLB", "TIGR01774"], "any"),
    ("甲酸脱氢 fdh/fdn", "dir1", ["FDH", "TIGR01591", "TIGR01553", "TIGR01582", "TIGR01583"], "any"),
    ("丁酸通路 buk", "dir1", ["TIGR02707"], "any"),
    ("丁酸通路 β-氧化", "dir1", ["PF00441", "PF02771", "PF00378", "PF16113",
                                "PF00725", "PF02737", "PF01012"], "majority"),
    ("丙酸通路 prpD", "dir1", ["TIGR02330"], "any"),
    ("丙酸通路 pccB/mmgA", "dir1", ["TIGR01117", "TIGR00641", "TIGR00642"], "any"),
    # ---------------- dir2 ----------------
    ("胆盐水解酶 BSH", "dir2", ["PF02275"], "any"),
    ("bai 操纵子（7α-脱羟）", "dir2", ["BAI_N", "BAI_H", "BAI_CD"], "majority"),  # core trio only
    ("7α-HSDH", "dir2", ["HSDH_7A"], "any"),
    ("7β-HSDH", "dir2", ["HSDH_7B"], "any"),
    ("12α-HSDH", "dir2", ["HSDH_12A"], "any"),
    # ---------------- dir3 ----------------
    ("色氨酸合成 trpA/trpB", "dir3", ["TIGR00262", "TIGR00263", "TIGR01415"], "majority"),
    ("色氨酸酶 TnaA（产吲哚）", "dir3", ["TIGR02617", "TIGR02618", "PF08053", "TIGR02616"], "any"),
    ("色氨酸转运 mtr/tnaB", "dir3", ["TIGR00837"], "any"),
    ("犬尿氨酸酶 KynU", "dir3", ["TIGR01814"], "any"),
    ("犬尿氨酸通路 KynA/B", "dir3", ["TIGR03035", "TIGR03036"], "any"),
    ("芳香族脱羧酶 TyrDC/Tdc", "dir3", ["TYDC", "TIGR03811"], "any"),
    ("谷氨酸脱羧酶 Gad", "dir3", ["TIGR03812"], "any"),
]

# load marker presence
rows = list(csv.DictReader(open(os.path.join(DM, "unit_family_matrix.tsv"),
                                encoding="utf-8"), delimiter="\t"))
# column header: "<PFxxxxx>|<hmm_name>" or "<family>|<hmm_name>"
col_of = {}
for c in rows[0]:
    if c.startswith("PF") or c.startswith("TIGR") or c in ("ADHE", "PFLB", "FDH", "LCTP",
                                                          "LLDEF", "TYDC", "ILD", "FLDC",
                                                          "BAI_N", "BAI_H", "BAI_CD", "BAI_F",
                                                          "BAI_K", "HSDH_7A", "HSDH_7B", "HSDH_12A"):
        col_of[c] = c

genus_of = {}
with open(os.path.join(ROOT, "metadata", "strain_manifest.tsv"), encoding="utf-8") as f:
    for r2 in csv.DictReader(f, delimiter="\t"):
        if r2.get("group") == "bacteria":
            genus_of[r2["unit_id"]] = (r2["current_name"].split()[0] if r2["current_name"] else "?")

missing = set()
for _, _, fams, _ in TRAITS:
    for f in fams:
        if f not in col_of:
            missing.add(f)
if missing:
    print("WARNING: markers absent from presence table:", sorted(missing))

with open(os.path.join(DM, "unit_functional_guilds.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["unit_id", "cn_name", "genus"] + [t[0] for t in TRAITS])
    for r in rows:
        vals = []
        for name, d, fams, rule in TRAITS:
            fs = [float(r.get(col_of[f], 0)) for f in fams if f in col_of]
            if not fs:
                vals.append(0.0)
                continue
            if rule == "any":
                v = max(fs)
            else:
                v = sum(1 for x in fs if x >= 0.5) / len(fs)
            vals.append(round(v, 2))
        w.writerow([r["unit_id"], r["cn_name"], genus_of.get(r["unit_id"], "?")] + vals)

# summary
lines = ["功能性状矩阵（阈值：单元内 ≥60%% 基因组携带记为阳性；%d 个单元）" % len(rows), ""]
by_trait = {t[0]: [] for t in TRAITS}
for r in rows:
    for name, d, fams, rule in TRAITS:
        fs = [float(r.get(col_of[f], 0)) for f in fams if f in col_of]
        if not fs:
            v = 0
        elif rule == "any":
            v = max(fs)
        else:
            v = sum(1 for x in fs if x >= 0.5) / len(fs)
        if v >= 0.6:
            by_trait[name].append(r["cn_name"])
for name in by_trait:
    hits = by_trait[name]
    lines.append("%-24s %2d/%d 单元: %s" % (name, len(hits), len(rows), "、".join(hits) if hits else "—"))

txt = "\n".join(lines)
with open(os.path.join(DM, "unit_guild_summary.txt"), "w", encoding="utf-8") as f:
    f.write(txt + "\n")
print(txt)
