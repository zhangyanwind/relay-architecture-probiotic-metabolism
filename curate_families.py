# -*- coding: utf-8 -*-
"""Curate the 97 direction HMM families into three interpretability classes.

  marker  the HMM pinpoints one pathway step (BSH, bai*, trpA/B, tnaA, buk,
          pta, ackA, larA, pflB, adhE, caiA/bcd, prpD, HSDH_7A ...)
  family  a pathway enzyme family: right pathway context, but several paralogs
          (LDH/MDH, ECH, ETF, thiolase ...)
  broad   a generic domain family (SDR, aminotransferase, ThDP decarboxylase,
          FGGY kinase ...) -- presence does NOT equal the named gene

Marker-gene share = fraction of named hits whose gene symbol is in the
expected marker set (evidence that the hits are the intended gene).
"""
import csv
import os
import re
from collections import Counter

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
ANN = os.path.join(ROOT, "annotations")
PG = os.path.join(ROOT, "pangenome")
DM = os.path.join(PG, "direction_matrices")

# family -> (class, expected marker gene regex)
CURATION = {
    # ---------------- dir1 ----------------
    "Pfam:PF00056": ("family", r"^(ldh|lctd|lcte|mdh)"),
    "Pfam:PF00378": ("family", r"^(crt|caid|ldef|fadb|menb)"),
    "Pfam:PF00441": ("family", r"^(caia|bcd|fadj|aco)"),
    "Pfam:PF00725": ("family", r"^(clbd|mmgb|fadb)"),
    "Pfam:PF01012": ("family", r"^(acra|acrb|fixa|etfb)"),
    "Pfam:PF01565": ("broad", r"^(murb)"),
    "Pfam:PF02615": ("broad", r"^(alld|ldha|ldh)"),   # hits allantoin DH, not D-LDH
    "Pfam:PF02652": ("marker", r"^(lutp|lldp|lctp)"),
    "Pfam:PF02737": ("family", r"^(clbd|mmgb|faib)"),
    "Pfam:PF02771": ("family", r"^(caia|bcd)"),
    "Pfam:PF02782": ("broad", r"^(glpk|gntk|xylb|rhb)"),
    "Pfam:PF02826": ("broad", r"^(sera|ghrb|ghra)"),
    "Pfam:PF09330": ("marker", r"^(dld)"),
    "Pfam:PF09861": ("marker", r"^(lara)"),
    "Pfam:PF16113": ("family", r"^(menb|caid|crt)"),
    "Pfam:PF21113": ("marker", r"^(lara)"),
    "TIGRFAMs:TIGR00016": ("marker", r"^(acka|tdcd)"),
    "TIGRFAMs:TIGR00641": ("family", r"^(scpa|mcm|mut)"),
    "TIGRFAMs:TIGR00642": ("family", r"^(scpa|muta|mcm)"),
    "TIGRFAMs:TIGR00651": ("marker", r"^(pta|ptb)"),
    "TIGRFAMs:TIGR00790": ("broad", r"^(nirc|foca)"),
    "TIGRFAMs:TIGR00795": ("marker", r"^(lutp|lldp|lctp)"),
    "TIGRFAMs:TIGR01108": ("family", r"^(pyc)"),
    "TIGRFAMs:TIGR01109": ("family", r"."),
    "TIGRFAMs:TIGR01117": ("marker", r"^(ldee|mmga|pccb)"),
    "TIGRFAMs:TIGR01235": ("family", r"^(pyc)"),
    "TIGRFAMs:TIGR01255": ("marker", r"^(pflb|dhab)"),
    "TIGRFAMs:TIGR01553": ("marker", r"^(fdng|fdog)"),
    "TIGRFAMs:TIGR01582": ("family", r"^(fdnh|fdoh|hyba)"),
    "TIGRFAMs:TIGR01583": ("marker", r"^(fdni|fdoi)"),
    "TIGRFAMs:TIGR01591": ("marker", r"^(fdhf|fdha|nas)"),
    "TIGRFAMs:TIGR01722": ("family", r"^(alda|adhe|prua)"),
    "TIGRFAMs:TIGR01756": ("family", r"^(mdh|ldh|lcte)"),
    "TIGRFAMs:TIGR01763": ("family", r"^(mdh|ldh|lcte)"),
    "TIGRFAMs:TIGR01771": ("family", r"^(mdh|ldh|lcte)"),
    "TIGRFAMs:TIGR01774": ("marker", r"^(pflb)"),
    "TIGRFAMs:TIGR01829": ("broad", r"^(fabg|fabi|ymfi)"),
    "TIGRFAMs:TIGR01930": ("family", r"^(fada|atob|fadi)"),
    "TIGRFAMs:TIGR02188": ("broad", r"^(dlta|acsa)"),
    "TIGRFAMs:TIGR02316": ("broad", r"^(dlta|acsa)"),
    "TIGRFAMs:TIGR02330": ("marker", r"^(prpd)"),
    "TIGRFAMs:TIGR02445": ("family", r"^(fada|atob)"),
    "TIGRFAMs:TIGR02446": ("family", r"^(fada|atob)"),
    "TIGRFAMs:TIGR02493": ("broad", r"^(nrdg|pfla)"),
    "TIGRFAMs:TIGR02518": ("marker", r"^(adhe|alda)"),
    "TIGRFAMs:TIGR02706": ("marker", r"^(pta|ptb)"),
    "TIGRFAMs:TIGR02707": ("marker", r"^(buk)"),
    "TIGRFAMs:TIGR02708": ("family", r"."),
    "TIGRFAMs:TIGR02720": ("broad", r"^(spxb|alss|ilvb|oxc)"),
    "TIGRFAMs:TIGR02751": ("family", r"^(ppca)"),
    "TIGRFAMs:TIGR03081": ("broad", r"^(gloa|mce)"),
    "TIGRFAMs:TIGR03215": ("family", r"^(mhpf)"),
    "TIGRFAMs:TIGR03458": ("family", r"^(scpc)"),
    "TIGRFAMs:TIGR03948": ("family", r"^(scpc)"),
    "TIGRFAMs:TIGR04060": ("broad", r"^(nirc|foca)"),
    "custom:ADHE": ("marker", r"^(adhe)"),
    "custom:FDH": ("marker", r"^(fdhf|fdha)"),
    "custom:LCTP": ("marker", r"^(lutp|lldp|lctp)"),
    "custom:LLDEF": ("family", r"."),
    "custom:PFLB": ("marker", r"^(pflb)"),
    # ---------------- dir2 ----------------
    "Pfam:PF01073": ("broad", r"^(gale|rfbb|rfbd|hsd)"),
    "Pfam:PF01758": ("broad", r"^(arsb)"),
    "Pfam:PF02275": ("marker", r"^(bsh|cbh|yxei)"),
    "Pfam:PF03446": ("broad", r"^(gnda|gnd)"),
    "TIGRFAMs:TIGR00841": ("family", r"^(arsb|yocs|bass)"),
    "custom:BAI_CD": ("marker", r"^(bilr|baicd|baib)"),
    "custom:BAI_F": ("marker", r"^(caib|baif)"),
    "custom:BAI_H": ("marker", r"^(bilr|baih)"),
    "custom:BAI_K": ("marker", r"^(caib|baik)"),
    "custom:BAI_N": ("marker", r"^(ytfp|bain)"),
    "custom:HSDH_12A": ("marker", r"^(hsd|hdh)"),
    "custom:HSDH_7A": ("marker", r"^(hdha|hsd)"),
    "custom:HSDH_7B": ("marker", r"^(hdhb|hsd)"),
    # ---------------- dir3 ----------------
    "Pfam:PF00155": ("broad", r"^(hisc|glya|pdxr|dapc)"),
    "Pfam:PF00221": ("family", r"^(huth|tna)"),
    "Pfam:PF01488": ("family", r"^(hema|aroe)"),
    "Pfam:PF08053": ("marker", r"^(tnac)"),
    "TIGRFAMs:TIGR00262": ("marker", r"^(trpa)"),
    "TIGRFAMs:TIGR00263": ("marker", r"^(trpb|ilva)"),
    "TIGRFAMs:TIGR00837": ("marker", r"^(mtr|tnab|tyrp)"),
    "TIGRFAMs:TIGR01415": ("marker", r"^(trpb)"),
    "TIGRFAMs:TIGR01814": ("marker", r"^(kynu)"),
    "TIGRFAMs:TIGR02616": ("marker", r"^(tnac)"),
    "TIGRFAMs:TIGR02617": ("marker", r"^(tnaa)"),
    "TIGRFAMs:TIGR02618": ("marker", r"^(tnaa)"),
    "TIGRFAMs:TIGR03035": ("marker", r"^(kynb)"),
    "TIGRFAMs:TIGR03036": ("marker", r"^(kyna)"),
    "TIGRFAMs:TIGR03393": ("broad", r"^(spxb|alss|ilvb|ipdc|poxb)"),
    "TIGRFAMs:TIGR03811": ("marker", r"^(tdc)"),
    "TIGRFAMs:TIGR03812": ("family", r"^(gad)"),
    "custom:FLDC": ("family", r"^(yjim)"),
    "custom:ILD": ("family", r"."),
    "custom:TYDC": ("marker", r"^(tdc)"),
}

hits = list(csv.DictReader(open(os.path.join(ANN, "hmm_hits_bacteria.tsv"),
                               encoding="utf-8"), delimiter="\t"))
by_fam = {}
for r in hits:
    by_fam.setdefault("%s:%s" % (r["source"], r["hmm_acc"]), []).append(r)

rows = []
for fam, lst in by_fam.items():
    cls, pattern = CURATION.get(fam, ("broad", r"."))
    named = [r["gene"] for r in lst if r["gene"] and r["gene"].upper() != "NA"]
    match = [g for g in named if re.match(pattern, g.lower())]
    top = Counter(named).most_common(3)
    rows.append({
        "family": fam,
        "direction": lst[0]["direction"],
        "hmm_name": lst[0]["hmm_name"],
        "class": cls,
        "marker_genes": pattern.strip("^$").replace("(", "").replace(")", "").replace("|", "/"),
        "n_hits": len(lst),
        "n_named": len(named),
        "marker_share": round(len(match) / len(named), 3) if named else "",
        "top_genes": ",".join("%s(%d)" % (g, n) for g, n in top),
    })

out = os.path.join(DM, "family_curation.tsv")
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in sorted(rows, key=lambda x: (x["direction"], {"marker": 0, "family": 1, "broad": 2}[x["class"]], -x["n_hits"])):
        w.writerow(r)

print("saved ->", out)
for d in ("dir1", "dir2", "dir3"):
    sub = [r for r in rows if r["direction"] == d]
    c = Counter(r["class"] for r in sub)
    print("%s: %d families | marker %d, family %d, broad %d" % (d, len(sub), c["marker"], c["family"], c["broad"]))
