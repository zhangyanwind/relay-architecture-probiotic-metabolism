# -*- coding: utf-8 -*-
"""Build custom HMMs for pathway families missing from TIGRFAMs/Pfam.

For each target family: fetch reviewed UniProt sequences via the REST API,
align with FAMSA (Windows binary), build a profile HMM with pyhmmer, and
collect all HMMs into one file for hmmscan/hmmsearch.

Outputs (under D:\\Probiotics_Pangenome\\tools\\hmm\\custom\\):
  <family>.faa          raw fetched sequences
  <family>.aln.faa      FAMSA alignment
  <family>.hmm          pyhmmer HMM
  custom_all.hmm        concatenated HMM library
  build_log.tsv         family -> n_sequences, status
"""
import os
import sys
import csv
import json
import time
import urllib.request
import urllib.parse
import subprocess

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
HMM_DIR = os.path.join(ROOT, "tools", "hmm")
CUSTOM_DIR = os.path.join(HMM_DIR, "custom")
FAMSA = os.path.join(ROOT, "tools", "famsa", "famsa.exe")
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
UA = {"User-Agent": "probiotic-pangenome-stage2/1.0"}
MAX_SEQ = 80          # cap per family (keeps FAMSA fast)
MIN_SEQ = 1           # thin families (1-2 seqs) are flagged, not fatal

# family key -> (UniProt query, human description)
FAMILIES = {
    # ---- direction 1: lactate -> SCFA ----
    "LCTP":   ('gene:lctP AND reviewed:true AND taxonomy_id:2', "lactate permease LctP"),
    "DLD":    ('gene:dldA AND reviewed:true AND taxonomy_id:2', "D-lactate dehydrogenase DldA"),
    "LLDEF":  ('gene:lldE AND reviewed:true AND taxonomy_id:2', "L-lactate dehydrogenase (cytochrome) LldEF"),
    "PFLB":   ('gene:pflB AND reviewed:true AND taxonomy_id:2', "pyruvate formate-lyase PflB"),
    "POR":    ('(gene:porA OR gene:porG) AND reviewed:true AND taxonomy_id:2', "pyruvate:ferredoxin oxidoreductase"),
    "ADHE":   ('gene:adhE AND reviewed:true AND taxonomy_id:2', "bifunctional acetaldehyde-CoA/alcohol dehydrogenase"),
    # ---- direction 2: bile acids ----
    "HSDH_3A": ('(gene:baiA2 OR "3-alpha-hydroxysteroid dehydrogenase") AND reviewed:true AND taxonomy_id:2', "3alpha-hydroxysteroid dehydrogenase"),
    "HSDH_7A": ('(gene:hdhA OR "7-alpha-hydroxysteroid dehydrogenase") AND reviewed:true AND taxonomy_id:2', "7alpha-hydroxysteroid dehydrogenase"),
    "HSDH_12A": ('"12-alpha-hydroxysteroid dehydrogenase" AND reviewed:true AND taxonomy_id:2', "12alpha-hydroxysteroid dehydrogenase"),
    "HSDH_7B": ('"7-beta-hydroxysteroid dehydrogenase" AND reviewed:true AND taxonomy_id:2', "7beta-hydroxysteroid dehydrogenase"),
    "BAI_B":  ('gene:baiB AND reviewed:true AND taxonomy_id:2', "bile acid CoA ligase BaiB"),
    "BAI_CD": ('gene:baiCD AND reviewed:true AND taxonomy_id:2', "baiCD fused protein"),
    "BAI_E":  ('gene:baiE AND reviewed:true AND taxonomy_id:2', "baiE protein"),
    "BAI_F":  ('gene:baiF AND reviewed:true AND taxonomy_id:2', "baiF protein"),
    "BAI_G":  ('gene:baiG AND reviewed:true AND taxonomy_id:2', "baiG protein"),
    "BAI_H":  ('gene:baiH AND reviewed:true AND taxonomy_id:2', "baiH protein"),
    "BAI_I":  ('gene:baiI AND reviewed:true AND taxonomy_id:2', "baiI protein"),
    "BAI_J":  ('gene:baiJ AND reviewed:true AND taxonomy_id:2', "baiJ protein"),
    "BAI_K":  ('gene:baiK AND reviewed:true AND taxonomy_id:2', "baiK protein"),
    "BAI_N":  ('gene:baiN AND reviewed:true AND taxonomy_id:2', "baiN protein"),
    "BAI_O":  ('gene:baiO AND reviewed:true AND taxonomy_id:2', "baiO protein"),
    # ---- direction 3: tryptophan ----
    "ILD":    ('(gene:ild OR "indolelactate dehydrogenase") AND reviewed:true AND taxonomy_id:2', "indolelactate dehydrogenase"),
    "FLDC":   ('gene:fldC AND reviewed:true AND taxonomy_id:2', "phenyllactate dehydratase large subunit FldC"),
    "FLDB":   ('gene:fldB AND reviewed:true AND taxonomy_id:2', "phenyllactate dehydratase small subunit FldB"),
    "ASTD":   ('gene:astD AND reviewed:true AND taxonomy_id:2', "aromatic aldehyde dehydrogenase AstD"),
    "TYDC":   ('(gene:tyrDC OR gene:ddc) AND reviewed:true AND taxonomy_id:2', "aromatic L-amino acid decarboxylase"),
    "FDH":    ('gene:fdhF AND reviewed:true AND taxonomy_id:2', "formate dehydrogenase H FdhF"),
}


def http_get(url, timeout=120, retries=4):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError("GET failed %s (%s)" % (url, last))


def fetch_family(fam, query):
    """Fetch FASTA for a family query with progressive relaxation.

    Strict query first (reviewed + bacteria); if it yields nothing, retry
    without the reviewed constraint, then without the taxonomy constraint.
    Returns (records, mode).
    """
    variants = [query,
                query.replace(" AND reviewed:true", ""),
                query.replace(" AND reviewed:true", "").replace(" AND taxonomy_id:2", "")]
    for mode, q in zip(("strict", "unreviewed", "any_taxon"), variants):
        url = "%s?query=%s&format=fasta&size=%d" % (
            UNIPROT, urllib.parse.quote(q), MAX_SEQ)
        raw = http_get(url)
        text = raw.decode("utf-8", "replace")
        if text.startswith("Error") or "<html" in text[:200].lower():
            raise RuntimeError("UniProt error: %s" % text[:200])
        recs = []
        cur_id, buf = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if cur_id:
                    recs.append((cur_id, "".join(buf)))
                cur_id = line
                buf = []
            else:
                buf.append(line.strip())
        if cur_id:
            recs.append((cur_id, "".join(buf)))
        if recs:
            return recs, mode
    return [], "none"


def write_fasta(recs, path):
    with open(path, "w", encoding="utf-8") as f:
        for hdr, seq in recs:
            f.write(hdr + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def run_famsa(inp, outp):
    r = subprocess.run([FAMSA, "-t", "1", inp, outp],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("FAMSA failed: %s" % (r.stderr[-300:] or r.stdout[-300:]))


def build_hmm(aln_path, hmm_path, fam_name, desc):
    from pyhmmer.plan7 import Builder, Background
    from pyhmmer.easel import MSAFile, Alphabet
    with open(aln_path, "rb") as f:
        msas = list(MSAFile(f, alphabet=Alphabet.amino()))
    if not msas:
        raise RuntimeError("no MSA read from %s" % aln_path)
    tmsa = msas[0]
    # name/description must be set on the TextMSA before digitizing
    tmsa.name = fam_name.encode()
    tmsa.description = ("%s | custom HMM for pathway annotation" % desc).encode()
    msa = tmsa.digitize(Alphabet.amino())
    builder = Builder(Alphabet.amino())
    background = Background(Alphabet.amino())
    hmm = builder.build_msa(msa, background)[0]
    with open(hmm_path, "wb") as f:
        hmm.write(f)
    return hmm


def main():
    os.makedirs(CUSTOM_DIR, exist_ok=True)
    log = []
    built = []
    for fam, (query, desc) in FAMILIES.items():
        faa = os.path.join(CUSTOM_DIR, fam + ".faa")
        aln = os.path.join(CUSTOM_DIR, fam + ".aln.faa")
        hmm = os.path.join(CUSTOM_DIR, fam + ".hmm")
        try:
            recs, mode = fetch_family(fam, query)
            if len(recs) < MIN_SEQ:
                raise RuntimeError("no sequences found")
            write_fasta(recs, faa)
            run_famsa(faa, aln)
            build_hmm(aln, hmm, fam, desc)
            built.append(hmm)
            st = "OK" if len(recs) >= 5 else "OK_THIN"
            if mode != "strict":
                st += "+%s" % mode
            log.append({"family": fam, "n_seq": len(recs), "desc": desc,
                        "status": st, "query": query})
            print("[%s] %-9s %3d seqs (%s)  %s" % (st, fam, len(recs), mode, desc), flush=True)
        except Exception as e:
            log.append({"family": fam, "n_seq": 0, "desc": desc,
                        "status": "FAIL: %s" % e, "query": query})
            print("[FAIL] %-9s %r" % (fam, e), flush=True)
        time.sleep(0.4)

    # concatenate all HMMs
    all_hmm = os.path.join(CUSTOM_DIR, "custom_all.hmm")
    with open(all_hmm, "wb") as out:
        for hmm in built:
            with open(hmm, "rb") as f:
                out.write(f.read())
    with open(os.path.join(CUSTOM_DIR, "build_log.tsv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["family", "n_seq", "desc", "status", "query"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(log)
    n_ok = sum(1 for r in log if r["status"].startswith("OK"))
    print("=== built %d/%d families -> %s ===" % (n_ok, len(FAMILIES), all_hmm))
    if n_ok < len(FAMILIES):
        sys.exit(2)


if __name__ == "__main__":
    main()
