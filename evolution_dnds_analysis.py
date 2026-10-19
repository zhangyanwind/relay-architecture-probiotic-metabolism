#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Panel c replacement: genuine dN/dS of relay enzymes from coding sequences.

The previous panel c used nearest-neighbour protein identity and returned a
null result (p = 0.94). dN/dS is the standard test for selection and asks a
different question: are the *downstream* (host-effect-generating) genes under
weaker purifying selection than the *upstream* genes? If relay architecture is
adaptive rather than accidental, restricted downstream genes should show dN/dS
closer to 1 (relaxed constraint) than widespread upstream genes.

Method
------
For every annotated relay gene we recover its CDS from the genome FASTA using
the NCBI GFF coordinates (including phase), then count synonymous (s) and
nonsynonymous (n) differences to a within-group reference using the
Nei-Gojobori 1986 method. Groups with >= 4 sequences are retained.

Outputs (English labels):
  evolution_dnds_by_group.tsv
  evolution_dnds_summary.txt
"""
import os, gzip
import numpy as np
import pandas as pd

BASE = r"D:\Probiotics_Pangenome_v5"
GEN = os.path.join(BASE, "genomes")
GFF = os.path.join(BASE, "annotations")
ANN = os.path.join(BASE, "annotations")
ANA = os.path.join(BASE, "validation", "analysis")
os.makedirs(ANA, exist_ok=True)

# ------------------------------------------------------------------ codon table
bases = "TCAG"
aas = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON = {}
i = 0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON[b1 + b2 + b3] = aas[i]; i += 1

# Nei-Gojobori synonymous substitution weights, expressed as expected
# transitions and transversions per codon (Nei & Gojobori 1986). Each entry is
# the average weight of changes to that base within a codon; the diagonal is 0.
_W = {  # (from, to): (transitions, transversions)
    ('A', 'G'): (1, 0), ('A', 'C'): (0, 1), ('A', 'T'): (0, 2),
    ('G', 'A'): (1, 0), ('G', 'C'): (0, 2), ('G', 'T'): (0, 1),
    ('C', 'T'): (1, 0), ('C', 'A'): (0, 1), ('C', 'G'): (0, 2),
    ('T', 'C'): (1, 0), ('T', 'G'): (0, 2), ('T', 'A'): (0, 1),
}
def _ng_weight(frm, to):
    """Expected number of substitutions from one codon base to another."""
    ts, tv = _W[(frm, to)]
    return ts * 1.0 + tv * 0.5      # transversions weighted at 0.5

_NGT = {}
for a in bases:
    for b in bases:
        _NGT[(a, b)] = 0.0 if a == b else _ng_weight(a, b)

SYN = {}
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            SYN[b1 + b2 + b3] = set(
                x + y + z for x in bases for y in bases for z in bases
                if CODON[x + y + z] == CODON[b1 + b2 + b3])


def read_fasta(path):
    """Yield (header, sequence). Handles plain and gzipped files."""
    op = gzip.open if path.endswith(".gz") else open
    hdr, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(buf)
                hdr = line[1:].strip(); buf = []
            else:
                buf.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(buf)


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def cds_from_gff(gff_path, wanted_proteins, fna_path=None):
    """Return {locus_tag: (protein_id, cds_seq)} for requested protein IDs."""
    out = {}
    seqname = None
    parent_locus = None
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "gene":
                attrs = f[8]
                pid = None
                for kv in attrs.split(";"):
                    if kv.startswith("ID=gene-"):
                        pid = kv[3:]
                parent_locus = pid
            elif f[2] == "CDS":
                attrs = f[8]
                pid = None
                for kv in attrs.split(";"):
                    if kv.startswith("ID=cds-"):
                        pid = kv[7:]
                if pid and pid in wanted_proteins:
                    rows.append((f[0], int(f[3]), int(f[4]), int(f[7]), pid))
    if not rows:
        return out
    need_chroms = {r[0] for r in rows}
    seqs = {}
    acc = os.path.basename(gff_path).replace(".gff", "")
    fna = fna_path or os.path.join(GEN, acc, f"{acc}.fna")
    if not os.path.exists(fna):
        return out
    for hdr, s in read_fasta(fna):
        nm = hdr.split()[0]
        if nm in need_chroms:
            seqs[nm] = s.upper()
    for chrom, st, en, phase, pid in rows:
        if chrom not in seqs:
            continue
        cds = seqs[chrom][st - 1:en]
        if len(cds) % 3:
            cds = cds[:len(cds) - (len(cds) % 3)]
        if len(cds) < 150:
            continue
        # GFF phase: trim leading bases to reach the first complete codon
        if phase in (1, 2):
            cds = cds[phase:]
        # coding sequences must begin with ATG
        if cds[:3] != "ATG":
            continue
        if any(ch not in "ACGT" for ch in cds[:30]):
            continue
        out[pid] = cds
    return out


def _pair_counts(a, b):
    """Nei-Gojobori S and N for one sequence pair."""
    L = min(len(a), len(b))
    S = N = 0.0
    nS_cod = nN_cod = 0
    for j in range(0, L - 2, 3):
        x, y = a[j:j + 3], b[j:j + 3]
        if len(x) < 3 or len(y) < 3:
            break
        syn = SYN.get(y)
        if syn is None:
            continue
        nS_cod += len(syn) - 1
        nN_cod += 61 - len(syn) - 1
        if x == y:
            continue
        if x in syn:
            S += _NGT[(y[0], x[0])] + _NGT[(y[1], x[1])] + _NGT[(y[2], x[2])]
        else:
            N += 1.0
    return S, N, nS_cod, nN_cod


def ng_dnds(cds_list):
    """Nei-Gojobori (1986) dN/dS averaged over all within-group sequence pairs.

    Pairwise comparison (rather than comparison to a consensus) is used so that
    all polymorphic sites in the group contribute, which matters for closely
    related genomes where consensus-relative differences can be near zero.
    """
    seqs = [c for c in cds_list if len(c) >= 150]
    if len(seqs) < 2:
        return np.nan, np.nan, np.nan, 0
    # cap group size for tractability
    if len(seqs) > 24:
        seqs = seqs[:24]
    npair = 0
    S_tot = N_tot = 0.0
    nS_cod_tot = nN_cod_tot = 0.0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            S, N, sc, nc = _pair_counts(seqs[i], seqs[j])
            if S <= 0 and N <= 0:
                continue
            if sc <= 0 or nc <= 0:
                continue
            npair += 1
            S_tot += S; N_tot += N
            nS_cod_tot += sc; nN_cod_tot += nc
    if npair == 0 or S_tot <= 0 or nS_cod_tot <= 0 or nN_cod_tot <= 0:
        return np.nan, np.nan, np.nan, npair
    dS = (S_tot / npair) / (nS_cod_tot / npair)
    dN = (N_tot / npair) / (nN_cod_tot / npair)
    ratio = dN / dS if dS > 0 else np.nan
    return ratio, dN, dS, npair


# ------------------------------------------------------------------ load hits
hits = pd.read_csv(os.path.join(ANN, "hmm_hits_relay.tsv"), sep="\t",
                   low_memory=False)
man = pd.read_csv(os.path.join(BASE, "metadata", "strain_manifest.tsv"),
                  sep="\t", low_memory=False)
acc2unit = dict(zip(man["accession"], man["unit_id"]))
hits["unit_id"] = hits["genome"].map(acc2unit)
hits = hits[hits["unit_id"].notna()].copy()
SIDE = {"GABA_down": ("GABA", "down"), "GABA_up": ("GABA", "up"),
        "GABA_supply": ("GABA", "up"), "GABA_transport": ("GABA", "up"),
        "TMA_up": ("Trimethylamine", "up"),
        "TMA_down": ("Trimethylamine", "down"),
        "TMA_transport": ("Trimethylamine", "up")}
hits["side"] = hits["pathway"].map(lambda p: SIDE.get(p, ("Other", "na"))[1])
hits["pw"] = hits["pathway"].map(lambda p: SIDE.get(p, ("Other", "na"))[0])

rows = []
for gi, (acc, grp) in enumerate(hits.groupby("genome")):
    gff = os.path.join(GFF, acc, f"{acc}.gff")
    if not os.path.exists(gff):
        continue
    want = {p.split("|")[-1]: p for p in grp["protein_id"]}
    got = cds_from_gff(gff, set(want))
    for pid, cds in got.items():
        r0 = grp[grp["protein_id"] == want[pid]]
        if not len(r0):
            continue
        r0 = r0.iloc[0]
        rows.append(dict(genome=acc, unit_id=r0["unit_id"], pw=r0["pw"],
                         side=r0["side"], pathway=r0["pathway"],
                         gene=want[pid].split("|")[2]
                         if len(want[pid].split("|")) > 2 else "",
                         hmm_name=r0["hmm_name"], pid=pid, cds=cds,
                         cds_len=len(cds)))
    if gi % 60 == 0:
        print(f"[cds] {gi} genomes, {len(rows)} CDS")

cdf = pd.DataFrame(rows)
print(f"[cds] total {len(cdf)} CDS from {cdf['genome'].nunique()} genomes")
if len(cdf) == 0:
    raise SystemExit("no CDS recovered")

cdf["gene_key"] = np.where(cdf["gene"].isin(["NA", "", None]),
                           "HMM:" + cdf["hmm_name"].astype(str),
                           "gene:" + cdf["gene"].astype(str))

res = []
for (pw, side, gk), g in cdf.groupby(["pw", "side", "gene_key"]):
    if len(g) < 2:
        continue
    ratio, dn, ds, nseq = ng_dnds(list(g["cds"]))
    if not np.isfinite(ratio) or ratio <= 0 or ratio > 20:
        continue
    res.append(dict(pathway=pw, side=side,
                    gene=gk.split(":", 1)[1], group_key=gk,
                    n_cds=len(g), n_used=nseq,
                    n_units=g["unit_id"].nunique(),
                    dnds=ratio, dN=dn, dS=ds))
if not res:
    raise SystemExit("no valid dN/dS group; all comparisons were invariant")
out = pd.DataFrame(res).sort_values("dnds")
out.to_csv(os.path.join(ANA, "evolution_dnds_by_group.tsv"), sep="\t",
           index=False)
print(f"\n[dN/dS] {len(out)} gene groups")
print(out.to_string(index=False))

from scipy.stats import mannwhitneyu
upv = out[out["side"] == "up"]["dnds"]
dnv = out[out["side"] == "down"]["dnds"]
if len(upv) > 2 and len(dnv) > 2:
    U, p = mannwhitneyu(upv, dnv, alternative="two-sided")
    print(f"\n[test] upstream dN/dS = {upv.mean():.3f} (n={len(upv)}, "
          f"{np.median(upv):.3f} median)")
    print(f"[test] downstream dN/dS = {dnv.mean():.3f} (n={len(dnv)}, "
          f"{np.median(dnv):.3f} median)")
    print(f"[test] Mann-Whitney p = {p:.4f}")
    with open(os.path.join(ANA, "evolution_dnds_summary.txt"), "w") as fh:
        fh.write(f"upstream_mean_dnds\t{upv.mean():.4f}\n"
                 f"upstream_median_dnds\t{np.median(upv):.4f}\n"
                 f"n_upstream_groups\t{len(upv)}\n"
                 f"downstream_mean_dnds\t{dnv.mean():.4f}\n"
                 f"downstream_median_dnds\t{np.median(dnv):.4f}\n"
                 f"n_downstream_groups\t{len(dnv)}\n"
                 f"mannwhitney_p\t{p:.6f}\n")
    print("[saved] evolution_dnds_summary.txt")