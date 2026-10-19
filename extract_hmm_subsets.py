# -*- coding: utf-8 -*-
"""Extract curated HMM subsets from TIGRFAMs_all.hmm and Pfam-A.hmm.

Streams each source file block-by-block (blocks end with '//'), keeps blocks
whose ACC line is in the curated list, and writes combined subset files:
  tools/hmm/tigr_subset.hmm   (~57 families, directions 1-3)
  tools/hmm/pfam_subset2.hmm  (additions beyond the original 15)
Also writes a manifest TSV: acc -> source, name, desc, direction.
"""
import os
import csv
import sys

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
HMM_DIR = os.path.join(ROOT, "tools", "hmm")

# ---------------- curated lists ----------------
TIGR = {
    # ---- direction 1: lactate -> SCFA ----
    "TIGR00016": "dir1", "TIGR00651": "dir1", "TIGR02188": "dir1",
    "TIGR01255": "dir1", "TIGR02493": "dir1", "TIGR01774": "dir1",
    "TIGR02720": "dir1", "TIGR02518": "dir1", "TIGR03215": "dir1",
    "TIGR01930": "dir1", "TIGR02445": "dir1", "TIGR02446": "dir1",
    "TIGR01829": "dir1", "TIGR02706": "dir1", "TIGR02707": "dir1",
    "TIGR03948": "dir1", "TIGR01756": "dir1", "TIGR01771": "dir1",
    "TIGR00795": "dir1", "TIGR02708": "dir1", "TIGR01235": "dir1",
    "TIGR02316": "dir1", "TIGR02330": "dir1", "TIGR00640": "dir1",
    "TIGR00641": "dir1", "TIGR00642": "dir1", "TIGR03081": "dir1",
    "TIGR01117": "dir1", "TIGR01108": "dir1", "TIGR01109": "dir1",
    "TIGR00790": "dir1", "TIGR04060": "dir1", "TIGR01553": "dir1",
    "TIGR01582": "dir1", "TIGR01583": "dir1", "TIGR01591": "dir1",
    "TIGR01722": "dir1", "TIGR03458": "dir1", "TIGR01763": "dir1",
    "TIGR02751": "dir1",
    # ---- direction 2: bile acids ----
    "TIGR00841": "dir2", "TIGR00805": "dir2",
    # ---- direction 3: tryptophan ----
    "TIGR02617": "dir3", "TIGR02616": "dir3", "TIGR00262": "dir3",
    "TIGR00263": "dir3", "TIGR03811": "dir3", "TIGR03812": "dir3",
    "TIGR02618": "dir3", "TIGR03393": "dir3", "TIGR01814": "dir3",
    "TIGR03035": "dir3", "TIGR03036": "dir3", "TIGR03037": "dir3",
    "TIGR00837": "dir3", "TIGR01415": "dir3", "TIGR03429": "dir3",
    "TIGR03430": "dir3",
}

PFAM_ADD = {
    # ---- direction 1 additions (butyrate route enzymes missing from TIGRFAMs) ----
    "PF00378": "dir1", "PF16113": "dir1",          # crotonase / enoyl-CoA hydratase (crt)
    "PF02737": "dir1", "PF00725": "dir1", "PF18321": "dir1",  # 3-hydroxyacyl-CoA dh (hbd)
    "PF00441": "dir1", "PF02770": "dir1", "PF02771": "dir1",  # acyl-CoA dh (bcd)
    "PF00766": "dir1", "PF01012": "dir1",          # electron transfer flavoprotein (etfAB)
    "PF02652": "dir1", "PF11870": "dir1",          # lactate permease / LutABC lactate utilisation
    "PF02615": "dir1",                             # malate/L-lactate dehydrogenase
    "PF09330": "dir1",                             # D-lactate dehydrogenase membrane-binding (dld)
    "PF09861": "dir1", "PF21113": "dir1",          # lactate racemase
    # ---- direction 2 additions ----
    "PF01758": "dir2",                             # sodium bile acid symporter (SBF)
    "PF04775": "dir2",                             # BAAT bile acid-CoA:amino acid N-acyltransferase
    # ---- direction 3 additions ----
    "PF03301": "dir3",                             # tryptophan 2,3-dioxygenase (kynA)
    "PF11991": "dir3",                             # tryptophan dimethylallyltransferase (dmaW)
    "PF08053": "dir3",                             # tryptophanase operon leader peptide (tnaC)
    "PF01231": "dir3",                             # indoleamine 2,3-dioxygenase (ido)
    "PF30897": "dir3", "PF31012": "dir3",          # indolepyruvate oxidoreductase (iorA/iorB)
}


def extract(src_path, wanted, out_path, source, manifest):
    """Stream HMM blocks; keep those whose ACC is in `wanted`."""
    kept = 0
    acc_re = None
    import re
    with open(src_path, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:
        block = []
        in_block = False
        for line in fin:
            if line.startswith("HMMER3") or (line.startswith("NAME") and not in_block):
                in_block = True
                block = [line]
                continue
            if not in_block:
                continue
            block.append(line)
            if line.startswith("//"):
                in_block = False
                text = "".join(block)
                m = re.search(r"^ACC\s+(\S+)", text, re.M)
                if m:
                    acc_full = m.group(1)
                    acc = acc_full.split(".")[0]
                    if acc in wanted:
                        fout.write(text)
                        kept += 1
                        nm = re.search(r"^NAME\s+(\S+)", text, re.M)
                        de = re.search(r"^DESC\s+(.*)$", text, re.M)
                        manifest.append({
                            "acc": acc, "source": source,
                            "name": nm.group(1) if nm else "",
                            "desc": (de.group(1) if de else "")[:100],
                            "direction": wanted[acc],
                        })
                block = []
    return kept


def main():
    manifest = []
    tigr_src = os.path.join(HMM_DIR, "TIGRFAMs_all.hmm")
    pfam_src = os.path.join(HMM_DIR, "Pfam-A.hmm")
    tigr_out = os.path.join(HMM_DIR, "tigr_subset.hmm")
    pfam_out = os.path.join(HMM_DIR, "pfam_subset2.hmm")

    n1 = extract(tigr_src, TIGR, tigr_out, "TIGRFAMs", manifest)
    print("TIGRFAMs: kept %d / %d wanted -> %s" % (n1, len(TIGR), tigr_out), flush=True)
    missing1 = set(TIGR) - {m["acc"] for m in manifest if m["source"] == "TIGRFAMs"}
    if missing1:
        print("  MISSING from TIGRFAMs_all.hmm: %s" % ", ".join(sorted(missing1)))

    n2 = extract(pfam_src, PFAM_ADD, pfam_out, "Pfam", manifest)
    print("Pfam: kept %d / %d wanted -> %s" % (n2, len(PFAM_ADD), pfam_out), flush=True)
    missing2 = set(PFAM_ADD) - {m["acc"] for m in manifest if m["source"] == "Pfam"}
    if missing2:
        print("  MISSING from Pfam-A.hmm: %s" % ", ".join(sorted(missing2)))

    man_path = os.path.join(HMM_DIR, "hmm_manifest.tsv")
    with open(man_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["acc", "source", "name", "desc", "direction"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(manifest)
    print("manifest: %d families -> %s" % (len(manifest), man_path))
    return 0 if (n1 == len(TIGR) and n2 == len(PFAM_ADD)) else 1


if __name__ == "__main__":
    sys.exit(main())
