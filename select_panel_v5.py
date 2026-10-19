# -*- coding: utf-8 -*-
"""Panel v5: diversity-first selection.

Allocation N per species stays span-driven (v4 table). WITHIN a species the
strains are chosen by ANI distance, not by metadata order:

  1. mandatory: type strain / RefSeq reference / representative
  2. greedy max-min: each further pick maximises the minimum ANI distance to
     already-selected strains (farthest-point sampling)
  3. near-clonal collapse: a candidate whose ANI to ANY selected strain is
     >= 99.9% is skipped (unless it is itself a mandatory type/reference)
  4. stops at N_target, or earlier if no distinct strains remain
     (flagged LOW_DIVERSITY)

Usage: python select_panel_v5.py [pool_tsv]
"""
import csv
import glob
import io
import os
import sys
import gzip
import statistics
from collections import defaultdict

import numpy as np

ROOT = "C:/Users/zhangyan/WorkBuddy/2026-10-17-17-38-14"
SEL = os.path.join(ROOT, "ncbi_selection")
POOL = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SEL, "pool_top12.tsv")
NEW_GEN = "D:/Probiotics_Pangenome/genomes_v3"
OLD_GEN = "D:/Probiotics_Pangenome/genomes"

K = 21
SKETCH = 1000
FB = np.full(256, 255, dtype=np.uint8)
for i, c in enumerate("ACGT"):
    FB[ord(c)] = i
    FB[ord(c.lower())] = i
BASE = np.uint64(4)
COMP = np.full(256, 255, dtype=np.uint8)
COMP[ord("A")] = 3; COMP[ord("C")] = 2; COMP[ord("G")] = 1; COMP[ord("T")] = 0
COMP[ord("a")] = 3; COMP[ord("c")] = 2; COMP[ord("g")] = 1; COMP[ord("t")] = 0


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    buf = []
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line[0] != ">":
                buf.append(line.strip())
    return "".join(buf)


def sketch_of(seq):
    a = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    codes = FB[a]
    n = len(codes) - K + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64)
    hf = np.zeros(n, dtype=np.uint64); ok = np.ones(n, dtype=bool)
    for i in range(K):
        hf = hf * BASE + codes[i:i + n].astype(np.uint64)
        ok &= codes[i:i + n] != 255
    rc = COMP[a][::-1]
    hr = np.zeros(n, dtype=np.uint64); okr = np.ones(n, dtype=bool)
    for i in range(K):
        hr = hr * BASE + rc[i:i + n].astype(np.uint64)
        okr &= rc[i:i + n] != 255
    h = np.minimum(hf, hr)[ok & okr]
    return np.sort(h)[:SKETCH] if len(h) > 0 else h


def ani(a, b):
    if len(a) == 0 or len(b) == 0:
        return 0.0
    shared = len(set(a.tolist()) & set(b.tolist()))
    if shared == 0:
        return 0.0
    j = shared / (2.0 * SKETCH - shared)
    return 1.0 + (1.0 / K) * np.log(2 * j / (1 + j))


def find_genome(acc):
    for root in (NEW_GEN, OLD_GEN):
        for ext in (".fna", ".fna.gz"):
            p = os.path.join(root, acc, acc + ext)
            if os.path.exists(p):
                return p
    return None


def num(x):
    try:
        return float(x)
    except Exception:
        return None


def norm(s):
    return "".join(c for c in (s or "").lower() if c.isalnum())


CACHE = {}
CACHE_PATH = os.path.join(SEL, "sketch_cache.npz")
NEW_CACHE = {}
if os.path.exists(CACHE_PATH):
    try:
        z = np.load(CACHE_PATH)
        for k in z.files:
            CACHE[k.replace("_", ".", 1) if False else k] = z[k]
        print("cache loaded:", len(CACHE))
    except Exception as e:
        print("cache load failed:", e)

rows = list(csv.DictReader(open(POOL, encoding="utf-8"), delimiter="\t"))
by = defaultdict(list)
for r in rows:
    by[r["taxid"]].append(r)

# span-based N (same table as v4)
def n_from_span(span):
    if span is None:
        return 5
    if span < 0.5:
        return 5
    if span < 1.0:
        return 6 + int(round((span - 0.5) / 0.5 * 2))
    if span < 2.0:
        return 8 + int(round((span - 1.0) / 1.0 * 4))
    if span < 3.5:
        return 12 + int(round((span - 2.0) / 1.5 * 6))
    return min(25, 18 + int(round((span - 3.5) / 1.0 * 3)))


def prio(r):
    return (r["type_material"] == "True", r["refseq_cat"] == "reference genome",
            r["refseq_cat"] == "representative genome",
            num(r["checkm_comp"]) or 0, -(num(r["checkm_cont"]) or 100), r["release"] or "")


panel, summary = [], []
for t, lst in sorted(by.items(), key=lambda kv: kv[0]):
    # dedup identical strain names, then load sketches
    best = {}
    for r in lst:
        k = norm(r["strain"]) or r.get("biosample") or r["acc"]
        if k not in best or prio(r) > prio(best[k]):
            best[k] = r
    cand = []
    for r in best.values():
        p = find_genome(r["gcf"] or r["acc"])
        if p:
            cand.append(r)
    if not cand:
        continue
    a = [num(r["ani"]) for r in cand if num(r["ani"]) is not None]
    span = (max(a) - min(a)) if len(a) >= 2 else None
    n_target = n_from_span(span)
    if a and min(a) < 96:
        n_target = min(n_target, 12)

    # sketches (cached)
    sk = {}
    for r in cand:
        acc = r["gcf"] or r["acc"]
        key = acc.replace(".", "_")
        if key in CACHE:
            sk[acc] = CACHE[key]
        elif False:
            pass
        else:
            sk[acc] = sketch_of(read_fasta(find_genome(acc)))
            CACHE[key] = sk[acc]
            NEW_CACHE[key] = sk[acc]

    # mandatory
    mand = [r for r in cand if r["type_material"] == "True" or
            r["refseq_cat"] in ("reference genome", "representative genome")]
    mand.sort(key=prio, reverse=True)
    selected = []
    for r in mand:
        if len(selected) >= n_target:
            break
        acc = r["gcf"] or r["acc"]
        # skip if near-clonal to an already selected mandatory
        if any(ani(sk[acc], sk[s["gcf"] or s["acc"]]) >= 0.999 for s in selected):
            continue
        selected.append(r)
    if not selected:
        selected.append(max(cand, key=prio))

    def near_clonal(r):
        acc = r["gcf"] or r["acc"]
        return max((ani(sk[acc], sk[s["gcf"] or s["acc"]]) for s in selected), default=0) >= 0.999

    # greedy max-min
    while len(selected) < n_target:
        pool = [r for r in cand if r not in selected and not near_clonal(r)]
        if not pool:
            break
        scored = []
        for r in pool:
            acc = r["gcf"] or r["acc"]
            dmin = min(ani(sk[acc], sk[s["gcf"] or s["acc"]]) for s in selected)
            scored.append((dmin, prio(r), r))
        scored.sort(key=lambda x: (x[0], x[1]))
        selected.append(scored[0][2])

    # pairwise ANI stats + 97% lineage count of the selection
    if len(selected) >= 2:
        accs = [r["gcf"] or r["acc"] for r in selected]
        pairs = [ani(sk[accs[i]], sk[accs[j]]) for i in range(len(accs)) for j in range(i + 1, len(accs))]
        min_pair = min(pairs)
        med_pair = statistics.median(pairs)
        n = len(accs)
        par = list(range(n))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for i in range(n):
            for j in range(i + 1, n):
                if ani(sk[accs[i]], sk[accs[j]]) >= 0.97:
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        par[ri] = rj
        n_lin = len(set(find(i) for i in range(n)))
    else:
        min_pair = med_pair = None
        n_lin = len(selected)
    flag = "LOW_DIVERSITY" if len(selected) < n_target else ""
    if len(selected) < 5:
        flag = (flag + ",SMALL_PANEL").strip(",")
    panel.extend(selected)
    summary.append((t, selected[0]["cn_name"], len(cand), round(span, 2) if span is not None else "-",
                    n_target, len(selected), n_lin,
                    round(min_pair, 3) if min_pair else "-",
                    round(med_pair, 3) if med_pair else "-", flag))

out = os.path.join(SEL, "panel_v5.tsv")
cols = ["taxid", "cn_name", "acc", "gcf", "strain", "isolation", "refseq_cat", "type_material",
        "checkm_comp", "checkm_cont", "n_contigs", "length", "n50", "gc", "ani", "release",
        "organism", "prov_note"]
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for r in sorted(panel, key=lambda x: (x["cn_name"], x["strain"])):
        w.writerow(r)

print("%-18s %4s %6s %4s %4s %4s %7s %7s %s" % ("物种", "候选", "跨度", "N", "实际", "谱系", "最远对", "中位对", "标记"))
tot = 0
for t, name, n_c, sp, n_t, n_s, n_l, mn, md, flag in sorted(summary, key=lambda x: -x[5]):
    tot += n_s
    print("%-18s %4d %6s %4d %4d %4d %7s %7s %s" % (name[:16], n_c, sp, n_t, n_s, n_l, mn, md, flag))
print("%-18s %4s %6s %4s %5d" % ("合计", "", "", "", tot))
if NEW_CACHE:
    old = {}
    if os.path.exists(CACHE_PATH):
        try:
            z = np.load(CACHE_PATH)
            old = {k: z[k] for k in z.files}
        except Exception:
            pass
    old.update(NEW_CACHE)
    np.savez_compressed(CACHE_PATH, **old)
    print("cache saved:", len(old))

print("\nsaved ->", out)
