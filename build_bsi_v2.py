#!/usr/bin/env python3
import csv
import gzip
import os
import urllib.request
from collections import defaultdict

DBROOT = os.path.expanduser("~/Downloads/bsi_database_v2")

SUMMARY = f"{DBROOT}/downloads/assembly_summary_refseq.txt"
NODES = f"{DBROOT}/taxonomy_full/nodes.dmp"
NAMES = f"{DBROOT}/taxonomy_full/names.dmp"

GENOMES = f"{DBROOT}/genomes"
MASTER = f"{DBROOT}/master_reference"
TAXOUT = f"{DBROOT}/taxonomy_wfcompatible"

FASTA = f"{MASTER}/bsi_panel_v2.fasta"
KRAKEN_FASTA = f"{MASTER}/bsi_panel_v2.kraken_headers.fasta"
REF2TAXID = f"{MASTER}/bsi_panel_v2.ref2taxid.tsv"
MANIFEST = f"{DBROOT}/selected_assemblies_v2.tsv"

TARGETS = [
    ("Staphylococcus aureus", 1280),
    ("Staphylococcus epidermidis", 1282),
    ("Staphylococcus haemolyticus", 1283),
    ("Staphylococcus hominis", 1290),
    ("Staphylococcus capitis", 29388),
    ("Staphylococcus lugdunensis", 28035),
    ("Escherichia coli", 562),
    ("Klebsiella pneumoniae", 573),
    ("Klebsiella oxytoca", 571),
    ("Enterobacter cloacae", 550),
    ("Citrobacter freundii", 546),
    ("Serratia marcescens", 615),
    ("Proteus mirabilis", 584),
    ("Morganella morganii", 582),
    ("Pseudomonas aeruginosa", 287),
    ("Acinetobacter baumannii", 470),
    ("Enterococcus faecalis", 1351),
    ("Enterococcus faecium", 1352),
    ("Streptococcus pneumoniae", 1313),
    ("Streptococcus pyogenes", 1314),
    ("Streptococcus agalactiae", 1311),
    ("Streptococcus dysgalactiae", 1334),
    ("Aggregatibacter aphrophilus", 732),
    ("Homo sapiens", 9606),
]

KEEP_RANKS = {
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
}

def parse_dmp(line):
    return [p.strip() for p in line.rstrip("\n").split("|")]

def read_nodes():
    d = {}
    with open(NODES, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip():
                continue
            p = parse_dmp(line)
            d[int(p[0])] = {"parent": int(p[1]), "rank": p[2]}
    return d

def read_names():
    names = defaultdict(list)
    sci = {}
    with open(NAMES, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip():
                continue
            p = parse_dmp(line)
            tid = int(p[0])
            names[tid].append((p[1], p[2], p[3]))
            if p[3] == "scientific name":
                sci[tid] = p[1]
    return names, sci

def lineage(tid, nodes):
    out = []
    seen = set()
    cur = tid
    while True:
        if cur in seen:
            raise RuntimeError(f"Taxonomy cycle at {tid}")
        seen.add(cur)
        out.append(cur)
        parent = nodes[cur]["parent"]
        if parent == cur:
            break
        cur = parent
    return list(reversed(out))

def load_summary():
    lines = []
    header = None
    with open(SUMMARY, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#assembly_accession"):
                header = line[1:]
            elif line.startswith("# assembly_accession"):
                header = line[2:]
            elif line.startswith("#"):
                continue
            else:
                lines.append(line)
    if header is None:
        raise RuntimeError("Could not parse assembly_summary header")
    return list(csv.DictReader([header] + lines, delimiter="\t"))

def choose_best(rows):
    def score(r):
        refcat = r.get("refseq_category", "").lower()
        level = r.get("assembly_level", "").lower()
        latest = 1 if r.get("version_status", "").lower() == "latest" else 0
        refscore = 5 if refcat == "reference genome" else 4 if refcat == "representative genome" else 1
        levelscore = {"complete genome": 5, "chromosome": 4, "scaffold": 3, "contig": 2}.get(level, 0)
        return (latest, refscore, levelscore, r.get("seq_rel_date", ""), r.get("assembly_accession", ""))
    rows = [r for r in rows if r.get("ftp_path", "") not in ("", "na")]
    rows.sort(key=score, reverse=True)
    return rows[0]

def download(url, dest):
    url = url.replace("ftp://", "https://")
    tmp = dest + ".tmp"
    print("  " + url)
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as out:
        while True:
            b = r.read(1024 * 1024)
            if not b:
                break
            out.write(b)
    os.replace(tmp, dest)

def append_genome(gz_path, taxid, fasta_out, kraken_out, map_out):
    n = 0
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="ignore") as inp, \
         open(fasta_out, "a") as fa, \
         open(kraken_out, "a") as kfa, \
         open(map_out, "a") as mp:
        header = None
        seq = []

        def flush(h, s):
            nonlocal n
            if h is None:
                return
            acc = h.split()[0]
            fa.write(f">{acc}\n")
            kfa.write(f">{acc}|kraken:taxid|{taxid}\n")
            for i in range(0, len(s), 80):
                chunk = s[i:i+80] + "\n"
                fa.write(chunk)
                kfa.write(chunk)
            mp.write(f"{acc}\t{taxid}\n")
            n += 1

        for line in inp:
            if line.startswith(">"):
                flush(header, "".join(seq))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        flush(header, "".join(seq))
    return n

def main():
    os.makedirs(GENOMES, exist_ok=True)
    os.makedirs(MASTER, exist_ok=True)
    os.makedirs(TAXOUT, exist_ok=True)

    for p in [FASTA, KRAKEN_FASTA, REF2TAXID, MANIFEST]:
        if os.path.exists(p):
            os.remove(p)

    nodes = read_nodes()
    names, sci = read_names()
    rows = load_summary()

    by_species_taxid = defaultdict(list)
    for r in rows:
        try:
            by_species_taxid[int(r["species_taxid"])].append(r)
        except Exception:
            pass

    selected = []
    for sp, tid in TARGETS:
        candidates = by_species_taxid.get(tid, [])
        if not candidates:
            raise RuntimeError(f"No RefSeq assembly found for {sp} taxid {tid}")
        selected.append((sp, tid, choose_best(candidates)))

    with open(MANIFEST, "w") as man:
        man.write("species_name\ttaxid\tassembly_accession\trefseq_category\tassembly_level\tseq_rel_date\tftp_path\n")
        for sp, tid, r in selected:
            ftp = r["ftp_path"].replace("ftp://", "https://")
            asm = ftp.rstrip("/").split("/")[-1]
            url = f"{ftp}/{asm}_genomic.fna.gz"
            dest = f"{GENOMES}/{asm}_genomic.fna.gz"
            print(f"\nDownloading {sp} ({tid})")
            download(url, dest)
            n = append_genome(dest, tid, FASTA, KRAKEN_FASTA, REF2TAXID)
            print(f"  imported sequences: {n}")
            man.write(f"{sp}\t{tid}\t{r['assembly_accession']}\t{r['refseq_category']}\t{r['assembly_level']}\t{r.get('seq_rel_date','')}\t{ftp}\n")

    keep_raw = set()
    for _, tid, _ in selected:
        keep_raw.update(lineage(tid, nodes))

    final_keep = {1}
    for tid in keep_raw:
        if tid == 1:
            final_keep.add(tid)
            continue
        rank = nodes[tid]["rank"]
        lin = lineage(tid, nodes)
        bacterial = 2 in lin
        if rank == "no rank":
            continue
        if bacterial and rank == "kingdom":
            continue
        if rank in KEEP_RANKS:
            final_keep.add(tid)

    patched_parent = {}
    for tid in final_keep:
        if tid == 1:
            patched_parent[tid] = 1
            continue
        cur_parent = nodes[tid]["parent"]
        while cur_parent not in final_keep:
            cur_parent = nodes[cur_parent]["parent"]
        patched_parent[tid] = cur_parent

    with open(f"{TAXOUT}/nodes.dmp", "w") as out:
        for tid in sorted(final_keep):
            parent = patched_parent[tid]
            rank = nodes[tid]["rank"]
            out.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")

    with open(f"{TAXOUT}/names.dmp", "w") as out:
        for tid in sorted(final_keep):
            for name, unique, cls in names[tid]:
                out.write(f"{tid}\t|\t{name}\t|\t{unique}\t|\t{cls}\t|\n")

    open(f"{TAXOUT}/merged.dmp", "w").close()
    open(f"{TAXOUT}/delnodes.dmp", "w").close()

    print("\nDONE master reference")
    print(FASTA)
    print(KRAKEN_FASTA)
    print(REF2TAXID)
    print(TAXOUT)
    print(MANIFEST)

if __name__ == "__main__":
    main()
