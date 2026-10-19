# -*- coding: utf-8 -*-
"""Stage-2 genome download for the probiotic pangenome project.

Reads D:\\Probiotics_Pangenome\\metadata\\strain_manifest.tsv and downloads
GENOME_FASTA + GENOME_GFF for every selected accession via the NCBI Datasets
API (batches of 25 accessions per request).

Output layout:
  D:\\Probiotics_Pangenome\\genomes\\<accession>\\<accession>.fna
  D:\\Probiotics_Pangenome\\annotations\\<accession>\\<accession>.gff
  D:\\Probiotics_Pangenome\\metadata\\download_log.tsv   (status per accession)

Features: md5 verification against the md5sum.txt shipped inside each zip,
resume support (completed accessions are skipped), retry with backoff.
"""
import urllib.request
import urllib.parse
import io
import zipfile
import os
import csv
import time
import sys
import hashlib

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
UA = {"User-Agent": "probiotic-pangenome-stage1/1.0"}
ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome")
MANIFEST = os.path.join(ROOT, "metadata", "strain_manifest.tsv")
GEN_DIR = os.path.join(ROOT, "genomes")
ANN_DIR = os.path.join(ROOT, "annotations")
LOG_TSV = os.path.join(ROOT, "metadata", "download_log.tsv")
BATCH = 25
RETRIES = 4
LOG = []


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.append(line)


def http_get(url, timeout=300, retries=RETRIES):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            log("  retry %d/%d after %r" % (i + 1, retries, e))
            time.sleep(3 * (i + 1))
    raise RuntimeError("GET failed %s (%s)" % (url, last))


def md5_of_bytes(b):
    return hashlib.md5(b).hexdigest()


def done_marker(acc):
    """Resume check: fna + gff exist and md5 verified marker present.

    Accepts OK_NO_GFF genomes (GenBank assemblies shipped without annotation):
    those have a .md5_ok marker containing 'no_gff' and no .gff file.
    """
    fna = os.path.join(GEN_DIR, acc, acc + ".fna")
    gff = os.path.join(ANN_DIR, acc, acc + ".gff")
    ok = os.path.join(GEN_DIR, acc, ".md5_ok")
    if os.path.exists(fna) and os.path.exists(ok):
        marker = open(ok, encoding="utf-8", errors="replace").read()
        if "no_gff" in marker:
            return True
        return os.path.exists(gff)
    return False


def process_zip(data, accs):
    """Extract fna/gff for the given accessions from a batch zip; verify md5."""
    z = zipfile.ZipFile(io.BytesIO(data))
    names = z.namelist()
    # parse md5sum.txt -- key by the FULL relative path listed in the file.
    # (basing on basename() is wrong: every genome's GFF is named genomic.gff,
    #  so basename keys collide and only the last entry survives.)
    expect = {}
    for n in names:
        if n.endswith("md5sum.txt"):
            for line in z.read(n).decode("utf-8", "replace").splitlines():
                parts = line.split()
                if len(parts) == 2:
                    expect[parts[1]] = parts[0].lower()
            break
    results = {}
    for acc in accs:
        fna_name = None
        gff_name = None
        for n in names:
            if "/%s/" % acc in n and n.endswith(".fna"):
                fna_name = n
            if "/%s/" % acc in n and n.endswith(".gff"):
                gff_name = n
        if not fna_name:
            results[acc] = ("MISSING_IN_ZIP", 0, 0)
            continue
        if not gff_name:
            # GenBank-source assemblies without NCBI annotation: keep the FASTA,
            # mark no-gff (downstream annotation falls back to self-annotation).
            fna_b = z.read(fna_name)
            if fna_name in expect and expect[fna_name] != md5_of_bytes(fna_b):
                results[acc] = ("MD5_MISMATCH", len(fna_b), 0)
                continue
            d1 = os.path.join(GEN_DIR, acc)
            os.makedirs(d1, exist_ok=True)
            with open(os.path.join(d1, acc + ".fna"), "wb") as f:
                f.write(fna_b)
            with open(os.path.join(d1, ".md5_ok"), "w") as f:
                f.write("fna %s\nno_gff\n" % md5_of_bytes(fna_b))
            results[acc] = ("OK_NO_GFF", len(fna_b), 0)
            continue
        fna_b = z.read(fna_name)
        gff_b = z.read(gff_name)
        fna_md5 = md5_of_bytes(fna_b)
        gff_md5 = md5_of_bytes(gff_b)
        ok = True
        if fna_name in expect and expect[fna_name] != fna_md5:
            ok = False
        if gff_name in expect and expect[gff_name] != gff_md5:
            ok = False
        if not ok:
            results[acc] = ("MD5_MISMATCH", len(fna_b), len(gff_b))
            continue
        d1 = os.path.join(GEN_DIR, acc)
        d2 = os.path.join(ANN_DIR, acc)
        os.makedirs(d1, exist_ok=True)
        os.makedirs(d2, exist_ok=True)
        with open(os.path.join(d1, acc + ".fna"), "wb") as f:
            f.write(fna_b)
        with open(os.path.join(d2, acc + ".gff"), "wb") as f:
            f.write(gff_b)
        with open(os.path.join(d1, ".md5_ok"), "w") as f:
            f.write("fna %s\ngff %s\n" % (fna_md5, gff_md5))
        results[acc] = ("OK", len(fna_b), len(gff_b))
    return results


def main():
    os.makedirs(GEN_DIR, exist_ok=True)
    os.makedirs(ANN_DIR, exist_ok=True)
    rows = []
    with open(MANIFEST, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("accession"):
                rows.append(r)
    # unique accessions (B. longum etc. should not repeat, but be safe)
    seen, accs = set(), []
    for r in rows:
        a = r["accession"].strip()
        if a and a not in seen:
            seen.add(a)
            accs.append(a)
    log("manifest: %d rows, %d unique accessions" % (len(rows), len(accs)))

    todo = [a for a in accs if not done_marker(a)]
    log("to download: %d (already done: %d)" % (len(todo), len(accs) - len(todo)))

    status = {}
    t0 = time.time()
    for i in range(0, len(todo), BATCH):
        batch = todo[i:i + BATCH]
        url = "%s/genome/accession/%s/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF" % (
            API, ",".join(urllib.parse.quote(a) for a in batch))
        log("batch %d/%d: %s ..." % (i // BATCH + 1, (len(todo) + BATCH - 1) // BATCH, batch[0]))
        try:
            data = http_get(url)
        except Exception as e:
            for a in batch:
                status[a] = ("BATCH_FAILED", 0, 0)
            log("  batch failed: %r" % e)
            continue
        try:
            res = process_zip(data, batch)
            status.update(res)
        except Exception as e:
            for a in batch:
                status[a] = ("ZIP_ERROR", 0, 0)
            log("  zip error: %r" % e)
        ok = sum(1 for a in batch if status.get(a, ("",))[0] == "OK")
        log("  -> OK %d / %d (elapsed %.0fs)" % (ok, len(batch), time.time() - t0))
        time.sleep(0.5)

    # write log for every accession in manifest
    with open(LOG_TSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "status", "fna_bytes", "gff_bytes"])
        for a in accs:
            st = status.get(a)
            if st is None:
                if done_marker(a):
                    marker = open(os.path.join(GEN_DIR, a, ".md5_ok"), encoding="utf-8",
                                  errors="replace").read()
                    st = ("OK_NO_GFF", 0, 0) if "no_gff" in marker else ("OK_cached", 0, 0)
                else:
                    st = ("NOT_ATTEMPTED", 0, 0)
            while len(st) < 3:
                st = tuple(st) + (0,)
            w.writerow([a, st[0], st[1], st[2]])

    n_ok = sum(1 for a in accs if done_marker(a))
    log("=== DONE: %d/%d accessions complete ===" % (n_ok, len(accs)))
    if n_ok < len(accs):
        bad = [a for a in accs if not done_marker(a)]
        log("incomplete: %s" % ", ".join(bad[:50]))
        sys.exit(2)


if __name__ == "__main__":
    main()
