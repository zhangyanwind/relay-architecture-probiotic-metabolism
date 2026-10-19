# -*- coding: utf-8 -*-
"""Native protein clustering for the probiotics pangenome.

Replaces vsearch/CD-HIT (GitHub unreachable from this network) with an
equivalent two-phase pipeline:

  Phase 1  exact-duplicate collapse (100% identity, same as CD-HIT -c 1.0)
  Phase 2  MinHash-style LSH prefilter:
             - k=5 amino-acid k-mers (kmers containing non-standard
               residues are skipped)
             - bottom-64 hash sketch per protein
             - 32 bands x 2 rows banding -> candidate pairs
  Phase 3  verification: global alignment (parasail Needleman-Wunsch,
             BLOSUM62, open=11/extend=1)
             identity = identical_positions / min(len_a, len_b)   (CD-HIT def.)
             coverage = alignment_length  / min(len_a, len_b)
             family edge  : identity >= 0.90 and coverage >= 0.80
             weak edge    : identity >= 0.50 and coverage >= 0.50
  Phase 4  connected components (single linkage) on the family-edge graph
             -> gene families; genome presence -> core/soft-core/shell/cloud

Usage:
    python cluster_proteins.py prefilter     # phase 1+2, write candidates.npy
    python cluster_proteins.py verify        # phase 3, write edges.tsv.gz
    python cluster_proteins.py families      # phase 4, write family tables
    python cluster_proteins.py all           # everything
"""
import os
import csv
import sys
import time
import struct
import argparse
from collections import defaultdict

import numpy as np

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
PROT_DIR = os.path.join(ROOT, "proteins")
MANIFEST = os.path.join(ROOT, "metadata", "strain_manifest.tsv")
OUT_DIR = os.path.join(ROOT, "pangenome")
CAND_FILE = os.path.join(OUT_DIR, "candidates.npy")
EDGE_FILE = os.path.join(OUT_DIR, "edges.tsv.gz")
ASSIGN_FILE = os.path.join(OUT_DIR, "gene_family_assignments.tsv")
FAMILY_FILE = os.path.join(OUT_DIR, "family_table.tsv")
REP_FASTA = os.path.join(OUT_DIR, "family_representatives.fasta")
UNIQ_FASTA = os.path.join(OUT_DIR, "unique_proteins.fasta")
UNIQ_INDEX = os.path.join(OUT_DIR, "unique_proteins_index.tsv")

KMER_K = 5
N_BAND = 32          # MinHash bands, independent permutations
BAND_ROWS = 2        # r=2 S-curve: recall ~100% at J=0.42 (~90% id)
SKETCH = N_BAND * BAND_ROWS
SALTS = np.random.default_rng(20261018).integers(
    0, 2 ** 63, size=N_BAND * BAND_ROWS, dtype=np.uint64)
BUCKET_CAP = 3000    # skip LSH buckets larger than this (low-complexity)
ID_FAMILY = 0.90     # family-level identity cutoff
COV_FAMILY = 0.80    # family-level coverage of shorter seq (CD-HIT -aS)
ID_WEAK = 0.50       # ortholog-level identity
COV_WEAK = 0.50
MIN_MATCHES = 40     # absolute identical residues; kills short-peptide/window artifacts
COV_MAX = 2.0        # alignment length / min(len); >2 means one seq is a window match

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)
CODEX = np.full(256, 255, dtype=np.uint8)
for _i, _c in enumerate(AA):
    CODEX[ord(_c)] = _i


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


# ---------------------------------------------------------------- phase 1
def load_and_dedup():
    """Load all bacterial .faa, collapse identical sequences.

    Returns (accs, pids, seqs) for unique proteins + member lists.
    """
    log("loading proteins from %s" % PROT_DIR)
    with open(MANIFEST, encoding="utf-8") as f:
        manifest = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
    bacteria = {a for a, r in manifest.items() if r.get("group") != "yeast"}
    faa_files = sorted(f for f in os.listdir(PROT_DIR) if f.endswith(".faa"))
    faa_files = [f for f in faa_files if f[:-4] in bacteria]
    log("bacterial genomes with faa: %d" % len(faa_files))

    seq2id = {}          # sequence -> unique index
    uniq_seqs = []
    members = []         # unique index -> list of (acc, pid)
    n_total = 0
    for fn in faa_files:
        acc = fn[:-4]
        pid = None
        buf = []
        with open(os.path.join(PROT_DIR, fn), "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    if pid is not None:
                        _add(buf, seq2id, uniq_seqs, members, acc, pid)
                        n_total += 1
                    hdr = line[1:].split()[0]
                    parts = hdr.split("|")
                    pid = parts[1] if len(parts) > 1 else parts[0]
                    buf = []
                else:
                    buf.append(line.strip())
            if pid is not None:
                _add(buf, seq2id, uniq_seqs, members, acc, pid)
                n_total += 1
    log("loaded %d proteins -> %d unique sequences" % (n_total, len(uniq_seqs)))
    return faa_files, uniq_seqs, members


def _add(buf, seq2id, uniq_seqs, members, acc, pid):
    seq = "".join(buf).upper()
    if len(seq) < 5:
        return
    idx = seq2id.get(seq)
    if idx is None:
        idx = len(uniq_seqs)
        seq2id[seq] = idx
        uniq_seqs.append(seq)
        members.append([])
    members[idx].append((acc, pid))


# ---------------------------------------------------------------- phase 2
def kmer_hashes(seq):
    """numpy vectorized rolling hash of k-mers (skips kmers w/ non-std aa)."""
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


def prefilter(uniq_seqs):
    log("phase 2: sketches + LSH banding (%d unique proteins)" % len(uniq_seqs))
    sketches = np.zeros((len(uniq_seqs), N_BAND * BAND_ROWS), dtype=np.uint64)
    t0 = time.time()
    for i, seq in enumerate(uniq_seqs):
        sketches[i] = minhash_bands(kmer_hashes(seq))
        if (i + 1) % 100000 == 0:
            log("  sketch %d/%d (%.0fs)" % (i + 1, len(uniq_seqs), time.time() - t0))
    # banding
    bands = defaultdict(list)
    for i in range(len(uniq_seqs)):
        sk = sketches[i]
        for b in range(N_BAND):
            vals = [int(sk[b * BAND_ROWS + r]) for r in range(BAND_ROWS)]
            key = struct.pack("<B" + "Q" * BAND_ROWS, b, *vals)
            bands[key].append(i)
    log("  %d band keys, max bucket %d" %
        (len(bands), max(len(v) for v in bands.values()) if bands else 0))
    # candidate pairs
    pairs = []
    skipped = 0
    for key, idxs in bands.items():
        if len(idxs) < 2:
            continue
        if len(idxs) > BUCKET_CAP:
            skipped += 1
            continue
        idxs = np.asarray(idxs, dtype=np.int64)
        # all pairs within bucket
        m = len(idxs)
        ii, jj = np.triu_indices(m, k=1)
        pairs.append(idxs[ii] * np.int64(len(uniq_seqs)) + idxs[jj])
    if pairs:
        cand = np.unique(np.concatenate(pairs))
    else:
        cand = np.zeros(0, dtype=np.int64)
    log("  candidates: %d pairs (skipped %d oversized buckets)" % (len(cand), skipped))
    np.save(CAND_FILE, cand)
    log("saved -> %s" % CAND_FILE)


# ---------------------------------------------------------------- phase 3
_SEQ_CACHE = {}


def load_unique_seqs():
    if "seqs" not in _SEQ_CACHE:
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
        _SEQ_CACHE["seqs"] = seqs
    return _SEQ_CACHE["seqs"]


def verify_chunk(args):
    lo, hi = args
    import parasail
    seqs = load_unique_seqs()
    out = []
    cand = np.load(CAND_FILE, mmap_mode="r")
    for p in cand[lo:hi]:
        i, j = divmod(int(p), len(seqs))
        a, b = seqs[i], seqs[j]
        if len(a) < 5 or len(b) < 5:
            continue
        r = parasail.nw_stats(a, b, 11, 1, parasail.blosum62)
        shorter = min(len(a), len(b))
        if r.length == 0:
            continue
        ident = r.matches / shorter
        cov = r.length / shorter
        if ident >= ID_WEAK and cov >= COV_WEAK and r.matches >= MIN_MATCHES and cov <= COV_MAX:
            out.append((i, j, round(ident, 4), round(cov, 4)))
    return out


def verify(nproc=12):
    import gzip
    import multiprocessing as mp
    log("phase 3: verifying candidates with parasail NW")
    cand = np.load(CAND_FILE)
    log("  %d candidate pairs" % len(cand))

    n = len(cand)
    chunk = max(1, (n + nproc - 1) // nproc)
    ranges = [(lo, min(lo + chunk, n)) for lo in range(0, n, chunk)]
    t0 = time.time()
    with mp.Pool(nproc) as pool:
        results = pool.map(verify_chunk, ranges)
    edges = [e for r in results for e in r]
    log("  verified in %.0fs -> %d weak-or-strong edges" % (time.time() - t0, len(edges)))
    with gzip.open(EDGE_FILE, "wt", encoding="utf-8") as f:
        f.write("i\tj\tidentity\tcoverage\n")
        for i, j, idt, cov in edges:
            f.write("%d\t%d\t%s\t%s\n" % (i, j, idt, cov))
    log("saved -> %s" % EDGE_FILE)


# ---------------------------------------------------------------- phase 4
def families():
    import gzip
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    log("phase 4: building gene families")
    # unique protein index: members per unique seq
    members = []
    with open(UNIQ_INDEX, encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            members.append(parts[2].split(",") if len(parts) > 2 and parts[2] else [])
    n = len(members)
    log("  %d unique proteins" % n)

    with gzip.open(EDGE_FILE, "rt", encoding="utf-8") as f:
        next(f)
        rows, cols, ids, covs = [], [], [], []
        for line in f:
            p = line.rstrip("\n").split("\t")
            i, j = int(p[0]), int(p[1])
            idt, cov = float(p[2]), float(p[3])
            rows.append(i); cols.append(j)
            ids.append(idt); covs.append(cov)
    log("  %d verified edges (id>=0.5, cov>=0.5)" % len(rows))

    # genome per unique protein = first member's acc
    genome_of = []
    for m in members:
        genome_of.append(m[0].split("|")[0] if m else "")
    genome_of = np.array(genome_of)

    def components(id_cut, cov_cut):
        sel = [k for k in range(len(rows)) if ids[k] >= id_cut and covs[k] >= cov_cut]
        data = np.ones(len(sel), dtype=np.int8)
        g = coo_matrix((data, ([rows[k] for k in sel], [cols[k] for k in sel])),
                       shape=(n, n), dtype=np.int8)
        return connected_components(g, directed=False)

    ncomp, labels = components(ID_FAMILY, COV_FAMILY)
    log("  family-level (id>=%.2f, cov>=%.2f): %d families" %
        (ID_FAMILY, COV_FAMILY, ncomp))

    # assignments: (acc, pid, family)
    with open(ASSIGN_FILE, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "protein_id", "family_id", "is_representative",
                    "seq_len", "n_members"])
        seq_lens = []
        with open(UNIQ_FASTA, encoding="utf-8") as ff:
            L = 0
            for line in ff:
                if line.startswith(">"):
                    if L:
                        seq_lens.append(L)
                    L = 0
                else:
                    L += len(line.strip())
            seq_lens.append(L)
        # representative = longest member of each family
        rep_of_family = {}
        for u in range(n):
            fam = labels[u]
            if fam not in rep_of_family or seq_lens[u] > seq_lens[rep_of_family[fam]]:
                rep_of_family[fam] = u
        rep_set = set(rep_of_family.values())
        for u in range(n):
            fam = labels[u]
            for acc, pid in (m.split("|") for m in members[u]):
                w.writerow([acc, pid, fam, 1 if u in rep_set else 0,
                            seq_lens[u], len(members[u])])

    # family table
    fam_genomes = defaultdict(set)
    fam_nprot = defaultdict(int)
    for u in range(n):
        fam = labels[u]
        fam_nprot[fam] += len(members[u])
        for m in members[u]:
            fam_genomes[fam].add(m.split("|")[0])
    with open(MANIFEST, encoding="utf-8") as f:
        manifest = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
    n_genomes_total = len(set(m.split("|")[0] for mem in members for m in mem))
    log("  total genomes in clustering: %d" % n_genomes_total)

    with open(FAMILY_FILE, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["family_id", "representative_protein", "n_proteins", "n_genomes",
                    "genome_fraction", "category", "genomes"])
        for fam in sorted(fam_nprot):
            g = fam_genomes[fam]
            frac = len(g) / n_genomes_total
            cat = ("core" if frac >= 0.99 else
                   "soft-core" if frac >= 0.95 else
                   "shell" if frac >= 0.15 else "cloud")
            rep_u = rep_of_family[fam]
            rep_acc, rep_pid = members[rep_u][0].split("|")
            w.writerow([fam, "%s|%s" % (rep_acc, rep_pid), fam_nprot[fam], len(g),
                        round(frac, 4), cat, ",".join(sorted(g))])
    log("saved -> %s, %s" % (ASSIGN_FILE, FAMILY_FILE))


def write_unique(uniq_seqs, members):
    with open(UNIQ_FASTA, "w", encoding="utf-8") as f, \
         open(UNIQ_INDEX, "w", encoding="utf-8") as fi:
        fi.write("uid\tn_members\tmembers\n")
        for i, (seq, mem) in enumerate(zip(uniq_seqs, members)):
            f.write(">u%d\n%s\n" % (i, seq))
            fi.write("%d\t%d\t%s\n" % (i, len(mem),
                                       ",".join("%s|%s" % m for m in mem)))


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prefilter", "verify", "families", "all"])
    ap.add_argument("--nproc", type=int, default=12)
    args = ap.parse_args()

    if args.stage in ("prefilter", "all"):
        if not os.path.exists(UNIQ_FASTA):
            faa_files, uniq_seqs, members = load_and_dedup()
            write_unique(uniq_seqs, members)
            log("saved -> %s, %s" % (UNIQ_FASTA, UNIQ_INDEX))
        else:
            log("reusing %s" % UNIQ_FASTA)
        # need uniq_seqs in memory for sketches
        uniq_seqs = []
        with open(UNIQ_FASTA, encoding="utf-8") as f:
            buf = []
            for line in f:
                if line.startswith(">"):
                    if buf:
                        uniq_seqs.append("".join(buf))
                    buf = []
                else:
                    buf.append(line.strip())
            if buf:
                uniq_seqs.append("".join(buf))
        prefilter(uniq_seqs)

    if args.stage in ("verify", "all"):
        verify(args.nproc)

    if args.stage in ("families", "all"):
        families()


if __name__ == "__main__":
    main()
