# -*- coding: utf-8 -*-
"""Second-pass orthogroup-level clustering (50% identity) via k=4 LSH.

Rationale: at the 90% identity level no family spans all 42 species
units (no protein-coding gene is >=90% identical across units), so the
collection-level core is empty and direction genes cannot be
distinguished from background. This pass uses k=4 k-mers (more sensitive
at moderate identity) to build orthogroups at >=50% identity / >=50%
coverage, the level at which conserved genes (ribosomal, translation)
do span units.

Stages:
  prefilter  k=4 sketches + LSH banding -> candidates_ortho.npy
  verify     verify NEW candidate pairs (skip pairs already in edges.tsv.gz)
             with the same crash-safe batched 3-worker design
  orthogroups  single-linkage components on the union edge set
  stats      orthogroup-level core/soft-core/shell/cloud + openness

Usage:
    python orthogroup_pass.py all
"""
import os
import csv
import time
import gzip
import struct
import argparse
from collections import defaultdict

import numpy as np

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
PG = os.path.join(ROOT, "pangenome")
CAND2 = os.path.join(PG, "candidates_ortho.npy")
UNIQ_FASTA = os.path.join(PG, "unique_proteins.fasta")
EDGES_90 = os.path.join(PG, "edges.tsv.gz")           # from the 90% pass
EDGES_NEW_RAW = os.path.join(PG, "edges_ortho_raw.tsv")
PROGRESS = os.path.join(PG, "verify_ortho_progress.txt")
EDGES_ALL = os.path.join(PG, "edges_all.tsv.gz")
ORTHO_ASSIGN = os.path.join(PG, "orthogroup_assignments.tsv")
ORTHO_TABLE = os.path.join(PG, "orthogroup_table.tsv")

KMER_K = 4
N_BAND = 32
BAND_ROWS = 2      # r=2: 3%/14%/45%/80% recall at J=0.03/0.07/0.14/0.22
SKETCH = N_BAND * BAND_ROWS
SALTS = np.random.default_rng(20261018).integers(
    0, 2 ** 63, size=N_BAND * BAND_ROWS, dtype=np.uint64)
BUCKET_CAP = 3000
BATCH = 200000
NPROC = int(os.environ.get("PROB_NPROC", "3"))
ID_CUT = 0.50
COV_CUT = 0.50
MIN_MATCHES = 40        # absolute identical residues; kills short-peptide/window artifacts
COV_MAX = 2.0           # alignment length / min(len); >2 means one seq is a window match

AA = "ACDEFGHIKLMNPQRSTVWY"
CODEX = np.full(256, 255, dtype=np.uint8)
for _i, _c in enumerate(AA):
    CODEX[ord(_c)] = _i

_SEQS = None


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def kmer_hashes(seq):
    b = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    codes = CODEX[b]
    n = len(codes) - KMER_K + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64)
    h = np.zeros(n, dtype=np.uint64)
    ok = np.ones(n, dtype=bool)
    for i in range(KMER_K):
        h = h * np.uint64(21) + codes[i:i + n].astype(np.uint64)
        ok &= codes[i:i + n] != 255
    return h[ok]


def load_seqs():
    global _SEQS
    if _SEQS is None:
        seqs = []
        with open(UNIQ_FASTA, encoding="utf-8") as f:
            buf = []
            for line in f:
                if line.startswith(">"):
                    if buf:
                        seqs.append("".join(buf))
                    buf = []
                else:
                    buf.append(line.strip())
            if buf:
                seqs.append("".join(buf))
        _SEQS = seqs
    return _SEQS


def mix64(x):
    """splitmix64 finalizer: strong avalanche so salted inputs give
    independent orderings (a mere multiplier perturbation does NOT)."""
    x = x ^ (x >> np.uint64(30))
    x = x * np.uint64(0xBF58476D1CE4E5B9)
    x = x ^ (x >> np.uint64(27))
    x = x * np.uint64(0x94D049BB133111EB)
    x = x ^ (x >> np.uint64(31))
    return x


def minhash_bands(h):
    """Standard MinHash LSH: one independent permutation per band.

    Band i value = min over k-mers of mix64(hash + salt_i). Two sequences
    share band i with P ~= J (k-mer Jaccard), independent across bands, so
    P(any of N_BAND bands) = 1-(1-J)^N_BAND.

    Two failure modes fixed here:
      * banding consecutive SORTED order statistics: one extra k-mer shifts
        every position -> all bands destroyed (47% recall bug);
      * multipliers differing by small constants (MIX_C + i): the orderings
        stay nearly identical, so all bands collapse into one.
    """
    out = np.empty(N_BAND * BAND_ROWS, dtype=np.uint64)
    if len(h) == 0:
        out[:] = np.iinfo(np.uint64).max
        return out
    for i in range(N_BAND * BAND_ROWS):
        out[i] = mix64(h + SALTS[i]).min()
    return out


# ------------------------------------------------------------ prefilter
def prefilter():
    seqs = load_seqs()
    n = len(seqs)
    log("k=%d prefilter over %d unique proteins" % (KMER_K, n))
    sketches = np.zeros((n, N_BAND * BAND_ROWS), dtype=np.uint64)
    t0 = time.time()
    for i, seq in enumerate(seqs):
        sketches[i] = minhash_bands(kmer_hashes(seq))
        if (i + 1) % 100000 == 0:
            log("  sketch %d/%d (%.0fs)" % (i + 1, n, time.time() - t0))
    bands = defaultdict(list)
    for i in range(n):
        sk = sketches[i]
        for b in range(N_BAND):
            vals = [int(sk[b * BAND_ROWS + r]) for r in range(BAND_ROWS)]
            key = struct.pack("<B" + "Q" * BAND_ROWS, b, *vals)
            bands[key].append(i)
    log("  %d band keys, max bucket %d" % (len(bands), max(len(v) for v in bands.values())))
    pairs = []
    skipped = 0
    for key, idxs in bands.items():
        if len(idxs) < 2:
            continue
        if len(idxs) > BUCKET_CAP:
            skipped += 1
            continue
        idxs = np.asarray(idxs, dtype=np.int64)
        m = len(idxs)
        ii, jj = np.triu_indices(m, k=1)
        pairs.append(idxs[ii] * np.int64(n) + idxs[jj])
    cand = np.unique(np.concatenate(pairs)) if pairs else np.zeros(0, dtype=np.int64)
    log("  candidates: %d pairs (skipped %d oversized buckets)" % (len(cand), skipped))
    np.save(CAND2, cand)

    # how many are already covered by the 90% pass?
    old = set()
    with gzip.open(EDGES_90, "rt", encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.split("\t")
            old.add(int(p[0]) * n + int(p[1]))
    new_mask = np.ones(len(cand), dtype=bool)
    # build set of new pairs via numpy searchsorted on sorted old keys
    old_arr = np.array(sorted(old), dtype=np.int64)
    idx = np.searchsorted(old_arr, cand)
    idx_clipped = np.clip(idx, 0, len(old_arr) - 1)
    new_mask = old_arr[idx_clipped] != cand
    log("  already verified in 90%% pass: %d, new to verify: %d" %
        ((~new_mask).sum(), new_mask.sum()))
    np.save(CAND2, cand[new_mask])
    log("saved -> %s (%d new pairs)" % (CAND2, len(cand[new_mask])))


# ------------------------------------------------------------ verify
def _init_worker():
    load_seqs()


def _verify_range(args):
    lo, hi = args
    import parasail
    seqs = _SEQS
    n = len(seqs)
    cand = np.load(CAND2, mmap_mode="r")
    out = []
    for p in cand[lo:hi]:
        i, j = divmod(int(p), n)
        a, b = seqs[i], seqs[j]
        la, lb = len(a), len(b)
        if la < 5 or lb < 5:
            continue
        r = parasail.nw_stats(a, b, 11, 1, parasail.blosum62)
        shorter = min(la, lb)
        if r.length == 0:
            continue
        ident = r.matches / shorter
        cov = r.length / shorter
        if ident >= ID_CUT and cov >= COV_CUT and r.matches >= MIN_MATCHES and cov <= COV_MAX:
            out.append("%d\t%d\t%.4f\t%.4f" % (i, j, ident, cov))
    return out


def verify():
    import multiprocessing as mp
    cand = np.load(CAND2)
    n = len(cand)
    done = 0
    if os.path.exists(PROGRESS):
        done = int(open(PROGRESS).read().strip() or 0)
    log("ortho verify: %d new pairs, done %d" % (n, done))
    if done >= n:
        return
    t0 = time.time()
    total = 0
    with mp.Pool(NPROC, initializer=_init_worker) as pool, \
         open(EDGES_NEW_RAW, "a", encoding="utf-8") as ef:
        pos = done
        bidx = 0
        while pos < n:
            hi = min(pos + BATCH, n)
            sub = max(1, (hi - pos + NPROC - 1) // NPROC)
            ranges = [(lo, min(lo + sub, hi)) for lo in range(pos, hi, sub)]
            results = pool.map(_verify_range, ranges)
            edges = [e for r in results for e in r]
            ef.write("\n".join(edges))
            if edges:
                ef.write("\n")
            ef.flush()
            os.fsync(ef.fileno())
            with open(PROGRESS, "w") as pf:
                pf.write(str(hi))
            pos = hi
            bidx += 1
            total += len(edges)
            el = time.time() - t0
            rate = (pos - done) / el if el > 0 else 0
            eta = (n - pos) / rate / 60 if rate > 0 else -1
            log("batch %d: %d/%d (%.0f%%), +%d edges, %.1f min, ETA %.1f min" %
                (bidx, pos, n, 100.0 * pos / n, len(edges), el / 60, eta))
    log("ortho verify done: %d new edges, %.1f min" % (total, (time.time() - t0) / 60))


# ------------------------------------------------------------ orthogroups
def orthogroups():
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    n = len(load_seqs())
    rows, cols = [], []
    with gzip.open(EDGES_90, "rt", encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.split("\t")
            rows.append(int(p[0])); cols.append(int(p[1]))
    n90 = len(rows)
    if os.path.exists(EDGES_NEW_RAW):
        with open(EDGES_NEW_RAW, encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue
                p = line.split("\t")
                rows.append(int(p[0])); cols.append(int(p[1]))
    log("union edges: %d (90%% pass %d + ortho pass %d)" % (len(rows), n90, len(rows) - n90))
    data = np.ones(len(rows), dtype=np.int8)
    g = coo_matrix((data, (rows, cols)), shape=(n, n), dtype=np.int8)
    ncomp, labels = connected_components(g, directed=False)
    log("orthogroups (id>=0.5, cov>=0.5, single linkage): %d" % ncomp)

    # members per unique protein
    members = []
    with open(os.path.join(PG, "unique_proteins_index.tsv"), encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            members.append(parts[2].split(",") if len(parts) > 2 and parts[2] else [])

    with open(ORTHO_ASSIGN, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "protein_id", "orthogroup_id"])
        for u in range(n):
            for m in members[u]:
                acc, pid = m.split("|")
                w.writerow([acc, pid, labels[u]])
    log("saved -> %s" % ORTHO_ASSIGN)

    fam_genomes = defaultdict(set)
    fam_nprot = defaultdict(int)
    for u in range(n):
        og = labels[u]
        fam_nprot[og] += len(members[u])
        for m in members[u]:
            fam_genomes[og].add(m.split("|")[0])
    with open(os.path.join(ROOT, "metadata", "strain_manifest.tsv"), encoding="utf-8") as f:
        manifest = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
    n_g = len(fam_genomes and set().union(*fam_genomes.values()))
    with open(ORTHO_TABLE, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["orthogroup_id", "n_proteins", "n_genomes", "genome_fraction",
                    "category", "genomes"])
        for og in sorted(fam_nprot):
            gset = fam_genomes[og]
            frac = len(gset) / n_g
            cat = ("core" if frac >= 0.99 else
                   "soft-core" if frac >= 0.95 else
                   "shell" if frac >= 0.15 else "cloud")
            w.writerow([og, fam_nprot[og], len(gset), round(frac, 4), cat,
                        ",".join(sorted(gset))])
    log("saved -> %s (%d orthogroups, %d genomes)" % (ORTHO_TABLE, len(fam_nprot), n_g))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prefilter", "verify", "orthogroups", "all"])
    args = ap.parse_args()
    if args.stage in ("prefilter", "all"):
        prefilter()
    if args.stage in ("verify", "all"):
        verify()
    if args.stage in ("orthogroups", "all"):
        orthogroups()


if __name__ == "__main__":
    main()
