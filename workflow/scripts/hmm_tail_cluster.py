"""Cluster HMMER tail-cassette hits per scaffold into candidate tailocin regions.

Reads:
  snakemake.input.tbl       : hmmsearch --tblout output (one row per gene-HMM hit)
  snakemake.input.gff       : Pyrodigal-gv GFF with gene coordinates
  snakemake.params.hmm_meta : detector-subset metadata (HMM name -> category)

Writes:
  snakemake.output.tsv : one row per candidate region with
                         scaffold, start, end, length, n_hits, n_distinct,
                         top_bitscore, hmms_present, category_set
"""
import csv, re, sys
from collections import defaultdict
from pathlib import Path

TBL_PATH  = Path(snakemake.input.tbl)
GFF_PATH  = Path(snakemake.input.gff)
HMM_META  = Path(snakemake.params.hmm_meta)
OUT_PATH  = Path(snakemake.output.tsv)
SAMPLE    = snakemake.params.sample
MAX_GAP_BP   = int(snakemake.params.max_gap_bp)
MIN_HITS     = int(snakemake.params.min_hits)
MIN_DISTINCT = int(snakemake.params.min_distinct)
MIN_BITS     = 30.0   # consistent with hmmsearch default + training-side filter

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# Detector-subset HMM names (cassette + lysis only, no capsid / portal)
detector_hmms = set()
hmm_cat = {}
with open(HMM_META) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        detector_hmms.add(r["hmm"])
        hmm_cat[r["hmm"]] = r.get("category", "")

# Gene coords from GFF
gene_coords = {}
with open(GFF_PATH) as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        attrs = dict(p.split("=", 1) for p in f[8].split(";") if "=" in p)
        gid = attrs.get("ID", "").rstrip(";")
        if not gid: continue
        try: s, e = int(f[3]), int(f[4])
        except ValueError: continue
        gene_coords[gid] = (f[0], min(s, e), max(s, e))

# Parse hmmsearch --tblout: per-gene top hit in detector subset
gene_top = {}   # gene_id -> (hmm_name, bitscore)
with open(TBL_PATH) as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.split()
        if len(f) < 6: continue
        target, _, query = f[0], f[1], f[2]
        try:
            score = float(f[5])
        except (IndexError, ValueError):
            continue
        if score < MIN_BITS: continue
        if query not in detector_hmms: continue
        prev = gene_top.get(target)
        if prev is None or score > prev[1]:
            gene_top[target] = (query, score)

# Cluster per-scaffold by gap-merging
by_scaffold = defaultdict(list)
for gid, (hmm, bits) in gene_top.items():
    coords = gene_coords.get(gid)
    if coords is None: continue
    sc, s, e = coords
    by_scaffold[sc].append((s, e, gid, hmm, bits))

regions = []
for sc, hits in by_scaffold.items():
    hits.sort(key=lambda x: x[0])
    cur = []
    last_end = None
    buckets = []
    for s, e, gid, hmm, bits in hits:
        if not cur:
            cur = [(s, e, gid, hmm, bits)]; last_end = e; continue
        if s - last_end <= MAX_GAP_BP:
            cur.append((s, e, gid, hmm, bits)); last_end = max(last_end, e)
        else:
            buckets.append(cur); cur = [(s, e, gid, hmm, bits)]; last_end = e
    if cur: buckets.append(cur)
    for c in buckets:
        distinct = {x[3] for x in c}
        if len(c) < MIN_HITS or len(distinct) < MIN_DISTINCT:
            continue
        rs = min(x[0] for x in c); re_ = max(x[1] for x in c)
        top_bits = max(x[4] for x in c)
        cats = sorted({hmm_cat.get(x[3], "") for x in c} - {""})
        regions.append((sc, rs, re_, re_ - rs, len(c), len(distinct),
                        round(top_bits, 2), ",".join(sorted(distinct)), ",".join(cats)))

with open(OUT_PATH, "w") as out:
    out.write("scaffold\tstart\tend\tlength\tn_hits\tn_distinct\ttop_bitscore\thmms_present\tcategory_set\n")
    for r in regions:
        out.write("\t".join(map(str, r)) + "\n")
print(f"[hmm_tail_cluster] {SAMPLE}: {len(regions)} candidate region(s)", file=sys.stderr)
