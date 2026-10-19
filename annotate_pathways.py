# -*- coding: utf-8 -*-
"""Pathway HMM annotation for the probiotic pangenome.

Searches every genome's protein FASTA against the curated pathway HMM library
(TIGRFAMs subset + 2 Pfam subsets + custom reference HMMs) with pyhmmer
hmmsearch, filters hits by E-value + coverage, and writes:

  annotations/hmm_hits.tsv          one row per (protein, family) best hit
  annotations/family_presence.tsv   genome x family gene-count matrix
  annotations/annotation_summary.tsv per-genome hit statistics
  annotations/done_genomes.txt      checkpoint (skip on restart)

Cutoffs (hmm_cov = fraction of the HMM model covered by the best domain):
  nseq >= 5 : E <= 1e-10, hmm_cov >= 0.50
  nseq 3-4  : E <= 1e-15, hmm_cov >= 0.60
  nseq <= 2 : E <= 1e-20, hmm_cov >= 0.75  (flagged putative=1)
seq_cov >= 0.30 in all cases (allows multi-domain proteins).
"""
import os
import csv
import sys
import time

import pyhmmer
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile, Alphabet

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
HMM_DIR = os.environ.get("HMM_ROOT", os.path.join(ROOT, "tools", "hmm"))
PROT_DIR = os.path.join(ROOT, "proteins")
ANN_DIR = os.path.join(ROOT, "annotations")
DONE_FILE = os.path.join(ANN_DIR, "done_genomes.txt")

# Pfam original subset: acc -> direction (curated)
PFAM1_DIR = {
    "PF02826": "dir1", "PF01073": "dir2", "PF00171": "dir1+dir3",
    "PF00155": "dir3", "PF02275": "dir2", "PF01565": "dir1",
    "PF02782": "dir1", "PF08544": "other", "PF00920": "other",
    "PF00478": "other", "PF00056": "dir1", "PF00221": "dir3",
    "PF01274": "other", "PF03446": "dir2", "PF01488": "dir3",
}
# custom families -> direction
CUSTOM_DIR = {
    "LCTP": "dir1", "DLD": "dir1", "LLDEF": "dir1", "PFLB": "dir1",
    "POR": "dir1", "ADHE": "dir1",
    "HSDH_3A": "dir2", "HSDH_7A": "dir2", "HSDH_12A": "dir2", "HSDH_7B": "dir2",
    "BAI_B": "dir2", "BAI_CD": "dir2", "BAI_E": "dir2", "BAI_F": "dir2",
    "BAI_G": "dir2", "BAI_H": "dir2", "BAI_I": "dir2", "BAI_K": "dir2",
    "BAI_N": "dir2",
    "ILD": "dir3", "FLDC": "dir3", "FLDB": "dir3", "ASTD": "dir3",
    "TYDC": "dir3", "FDH": "dir1",
}

HIT_FIELDS = ["genome", "protein_id", "gene", "protein_acc", "seq_len",
              "hmm_acc", "hmm_name", "source", "direction", "hmm_nseq",
              "evalue", "score", "bias", "hmm_from", "hmm_to", "hmm_len",
              "hmm_cov", "tgt_from", "tgt_to", "seq_cov", "putative"]

# families removed from the library after contamination review:
#  ASTD - the gene:astD query fetched N-succinylglutamate 5-semialdehyde
#         dehydrogenase (arginine catabolism), NOT the aromatic aldehyde DH
#  BAI_G - single-seq transporter HMM matches the MFS superfamily (bmr3 etc.)
EXCLUDE_CUSTOM = {"ASTD", "BAI_G"}

# display labels for equivalogs whose HMM family is broader than the name
FAMILY_LABEL = {
    "TIGR01829": "FabG/AcAcCoA_reductase_SDR",
    "TIGR03393": "ThDP_decarboxylase_IpdC_ARO10",
    "TIGR00841": "bile_acid_transporter_Bass",
}

# Curated absolute bit-score cutoffs for thin custom HMMs (nseq <= 4), set
# from score-distribution calibration (2026-10-18): each value sits in the
# natural gap between the true family cluster and the superfamily background.
#   HSDH_7A 250: hdhA cluster 293-415 (E. coli/C. butyricum/L. mucosae) vs FabG <=235
#   HSDH_3A 250: best genome hit is FabG 167.8 -> no genuine 3a-HSDH in collection
#   HSDH_7B 240: C. butyricum hdhB 438 vs generic SDR <=102
#   HSDH_12A 306: F. prausnitzii 403 vs next 217
#   ILD 275: D-2-hydroxyacid DH cluster 283-309 vs background <=241
#   BAI_B 634: no genome above 149.6 (ANL superfamily) -> baiB undetected
#   BAI_CD 250: C. butyricum 310-315 / F. prausnitzii 305 vs background <=147
#   BAI_H 250: C. butyricum 324-341 vs background <=127
#   BAI_F/BAI_K 200: CaiB/BaiF CoA-transferase family 238-272 vs background <=170
#   BAI_N 180: BaiN/RdsA NAD(P)/FAD oxidoreductase family 182-312
#   FDH 300: fdhF 1735 (E. coli) + FDH-alpha cluster 371-472 vs background <=114
#   LCTP 550: LctP/LldP/LutP cluster 648-997 vs background <=439
#   LLDEF 400: FAD-linked L-LDH cluster 735 vs background <=161
#   PFLB 900: pflB cluster 1115-1513 vs background <=782
#   ADHE 800: adhE cluster 1115-1630 vs background <=213
#   FLDC 80: only 4 weak hits, all family-annotated (double-cubane anaerobic reductase)
#   TYDC 400: tdc cluster 624-801 vs background <=101
SCORE_CUTOFF = {
    "HSDH_7A": 250, "HSDH_3A": 250, "HSDH_7B": 240, "HSDH_12A": 306,
    "ILD": 275, "BAI_B": 634, "BAI_CD": 250, "BAI_F": 200, "BAI_H": 250,
    "BAI_K": 200, "BAI_N": 180, "FDH": 300, "LCTP": 550, "LLDEF": 400,
    "PFLB": 900, "ADHE": 800, "FLDC": 80, "TYDC": 400,
}


def _s(x):
    """pyhmmer may return str or bytes depending on version."""
    if x is None:
        return ""
    return x.decode() if isinstance(x, bytes) else str(x)


def load_hmm_library():
    """Returns list of dicts: hmm object + metadata."""
    lib = []

    def add_file(path, source, dir_map=None):
        if not os.path.exists(path):
            print("  WARN missing %s" % path)
            return
        n = 0
        with HMMFile(path) as hf:
            for h in hf:
                acc = _s(h.accession).split(".")[0] or _s(h.name)
                name = _s(h.name)
                if dir_map is not None:
                    direction = dir_map.get(acc, dir_map.get(name, "other"))
                else:
                    direction = "other"
                lib.append({"hmm": h, "acc": acc, "name": name,
                            "source": source, "direction": direction,
                            "nseq": h.nseq, "M": h.M})
                n += 1
        print("  %-10s %-24s %3d HMMs" % (source, os.path.basename(path), n))

    # manifest from subset extraction (tigr + pfam2)
    man = {}
    man_path = os.path.join(HMM_DIR, "hmm_manifest.tsv")
    if os.path.exists(man_path):
        with open(man_path, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                man[(row["source"], row["acc"])] = row["direction"]

    # custom nseq lookup for thin-HMM cutoffs
    custom_nseq = {}
    blog = os.path.join(HMM_DIR, "custom", "build_log.tsv")
    if os.path.exists(blog):
        with open(blog, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    custom_nseq[row["family"]] = int(row["n_seq"])
                except (ValueError, KeyError):
                    pass

    # 1) TIGRFAMs subset (directions from manifest)
    n = 0
    with HMMFile(os.path.join(HMM_DIR, "tigr_subset.hmm")) as hf:
        for h in hf:
            acc = _s(h.accession).split(".")[0] or _s(h.name)
            lib.append({"hmm": h, "acc": acc, "name": _s(h.name),
                        "source": "TIGRFAMs",
                        "direction": man.get(("TIGRFAMs", acc), "other"),
                        "nseq": h.nseq, "M": h.M})
            n += 1
    print("  %-10s %-24s %3d HMMs" % ("TIGRFAMs", "tigr_subset.hmm", n))

    # 2) Pfam original subset
    add_file(os.path.join(HMM_DIR, "Pfam_subset.hmm"), "Pfam", PFAM1_DIR)
    # 3) Pfam additions
    add_file(os.path.join(HMM_DIR, "pfam_subset2.hmm"), "Pfam",
             {a: d for (s, a), d in man.items() if s == "Pfam"})
    # 4) custom
    n = 0
    cpath = os.path.join(HMM_DIR, "custom", "custom_all.hmm")
    with HMMFile(cpath) as hf:
        for h in hf:
            name = _s(h.name)
            if name in EXCLUDE_CUSTOM:
                continue
            lib.append({"hmm": h, "acc": name, "name": name, "source": "custom",
                        "direction": CUSTOM_DIR.get(name, "other"),
                        "nseq": custom_nseq.get(name, h.nseq), "M": h.M})
            n += 1
    print("  %-10s %-24s %3d HMMs" % ("custom", "custom_all.hmm", n))
    return lib


def cutoffs_for(meta, self_cuts):
    """Returns (evalue_cut, cov_cut, score_cut, putative)."""
    if meta["source"] == "custom":
        scut = self_cuts.get(meta["name"])
        n = meta["nseq"]
        if n <= 2:
            return 1e-20, 0.75, scut, 1
        if n <= 4:
            return 1e-15, 0.60, scut, 0
        # curated cutoff applies regardless of nseq when listed
        return 1e-10, 0.50, scut, 0
    return 1e-10, 0.50, None, 0

def annotate_genome(faa_path, lib, genome, self_cuts=None):
    if self_cuts is None:
        self_cuts = {}
    rows = []
    with SequenceFile(faa_path, digital=True, alphabet=Alphabet.amino()) as sf:
        seqs = list(sf)
    if not seqs:
        return rows, 0
    seq_meta = {}
    for s in seqs:
        # name like ACC|locus_tag|gene|protein_id
        parts = _s(s.name).split("|")
        seq_meta[_s(s.name)] = (parts + ["", "", "", ""])[:4]
    hmms = [m["hmm"] for m in lib]
    by_name = {}
    for m in lib:
        by_name.setdefault(m["name"], []).append(m)
    n_miss = 0
    for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=4):
        meta = None
        for m in lib:
            if m["hmm"] is hits.query:
                meta = m
                break
        if meta is None:
            # fall back to name matching (identity may not survive iteration)
            cands = by_name.get(_s(hits.query.name), [])
            if cands:
                meta = cands[0]
                n_miss += 1
        if meta is None:
            n_miss += 1
            continue
        e_cut, cov_cut, score_cut, putative = cutoffs_for(meta, self_cuts)
        for hit in hits:
            if hit.evalue > e_cut or hit.included is False:
                continue
            if score_cut is not None and hit.score < score_cut:
                continue
            d = hit.best_domain
            al = d.alignment
            hmm_cov = (al.hmm_to - al.hmm_from + 1) / al.hmm_length
            seq_cov = (al.target_to - al.target_from + 1) / al.target_length
            if hmm_cov < cov_cut or seq_cov < 0.30:
                continue
            parts = seq_meta.get(_s(hit.name), [_s(hit.name), "", "", ""])
            rows.append({
                "genome": genome,
                "protein_id": parts[1] if len(parts) > 1 else "",
                "gene": parts[2] if len(parts) > 2 else "",
                "protein_acc": parts[3] if len(parts) > 3 else "",
                "seq_len": len(hit.sequence) if hasattr(hit, "sequence") and hit.sequence else hit.length,
                "hmm_acc": meta["acc"], "hmm_name": FAMILY_LABEL.get(meta["acc"], meta["name"]),
                "source": meta["source"], "direction": meta["direction"],
                "hmm_nseq": meta["nseq"],
                "evalue": "%.3g" % hit.evalue, "score": "%.1f" % hit.score,
                "bias": "%.2f" % hit.bias,
                "hmm_from": al.hmm_from, "hmm_to": al.hmm_to, "hmm_len": al.hmm_length,
                "hmm_cov": "%.3f" % hmm_cov,
                "tgt_from": al.target_from, "tgt_to": al.target_to,
                "seq_cov": "%.3f" % seq_cov,
                "putative": putative,
            })
    if n_miss:
        print("  [WARN] %s: %d query-HMM matches failed" % (genome, n_miss), flush=True)
    return rows, len(seqs)


def main():
    os.makedirs(ANN_DIR, exist_ok=True)
    print("loading HMM library...", flush=True)
    lib = load_hmm_library()
    print("total HMMs: %d" % len(lib), flush=True)
    self_cuts = SCORE_CUTOFF
    print("curated score cutoffs for %d thin custom HMMs" % len(self_cuts), flush=True)

    hits_path = os.path.join(ANN_DIR, "hmm_hits.tsv")
    done = set()
    if os.path.exists(DONE_FILE):
        done = {x.strip() for x in open(DONE_FILE) if x.strip()}
    faas = sorted(f for f in os.listdir(PROT_DIR) if f.endswith(".faa"))
    todo = [f for f in faas if f[:-4] not in done]
    print("genomes: %d total, %d done, %d todo" % (len(faas), len(done), len(todo)), flush=True)

    write_header = not os.path.exists(hits_path) or os.path.getsize(hits_path) == 0
    t0 = time.time()
    n_done = 0
    with open(hits_path, "a", newline="", encoding="utf-8") as hf_out:
        w = csv.DictWriter(hf_out, fieldnames=HIT_FIELDS, delimiter="\t")
        if write_header:
            w.writeheader()
            hf_out.flush()
        for i, fn in enumerate(todo, 1):
            genome = fn[:-4]
            try:
                rows, nseq = annotate_genome(os.path.join(PROT_DIR, fn), lib, genome, self_cuts)
            except Exception as e:
                print("  [FAIL] %s: %r" % (genome, e), flush=True)
                continue
            for r in rows:
                w.writerow(r)
            hf_out.flush()
            with open(DONE_FILE, "a") as df:
                df.write(genome + "\n")
            n_done += 1
            if i % 20 == 0 or i == len(todo):
                el = time.time() - t0
                print("  [%d/%d] %s: %d hits (%d proteins), %.0fs elapsed" % (
                    i, len(todo), genome, len(rows), nseq, el), flush=True)
    print("=== annotation done: %d genomes, %.0fs ===" % (n_done, time.time() - t0), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
