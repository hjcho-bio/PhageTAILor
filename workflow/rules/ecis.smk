"""Module: ecis — cluster eCIStem MMseqs hits per scaffold into eCIS candidate regions.

Same per-scaffold gap-merging clusterer as `t6ss.smk`, but the reference DB is
the bundled eCIStem v1 protein set (`resources/ecistem_db/`). Mirrors the
training-side detector exactly so that user inference uses the same candidate
distribution the model was trained on.
"""

rule ecis_cluster_ecistem:
    benchmark:
        OUT/"benchmarks"/"ecis_cluster_ecistem.{s}.tsv"
    input:
        m8  = OUT/"annotation"/"{s}.ecistem.m8",
        gff = OUT/"annotation"/"{s}.gff",
    output:
        flag = touch(OUT/"ecis"/"{s}"/"_done"),
        tsv  = OUT/"ecis"/"{s}"/"candidates.tsv",
    conda: "../envs/python.yaml"
    params:
        meta = config["ecistem_metadata"],
        sample = lambda wc: wc.s,
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        python3 - <<PY
import csv, re
from collections import defaultdict

MIN_BITS, MIN_PIDENT = 50, 30
MAX_GAP_BP = 15000
MIN_DISTINCT, MIN_HITS = 3, 4

# Map eCIStem seqid to cluster_id; m8 target field has '|cluster=...' suffix
target_to_cluster = {{}}
with open("{params.meta}") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd: target_to_cluster[r["seqid"]] = r["cluster_id"]

def strip_ecistem_suffix(t): return t.split("|", 1)[0]

genes = {{}}
with open("{input.gff}") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        attrs = dict(p.split("=",1) for p in f[8].split(";") if "=" in p)
        gid = attrs.get("ID","").rstrip(";")
        if not gid: continue
        try: s, e = int(f[3]), int(f[4])
        except ValueError: continue
        genes[gid] = (f[0], min(s,e), max(s,e))

seen = set()
by_scaffold = defaultdict(list)
with open("{input.m8}") as fh:
    for line in fh:
        f = line.rstrip().split("\t")
        if len(f) < 6: continue
        q, t = f[0], f[1]
        try: pid, bits = float(f[2]), float(f[5])
        except ValueError: continue
        if bits < MIN_BITS or pid < MIN_PIDENT: continue
        if q in seen: continue
        seen.add(q)
        if q not in genes: continue
        cl = target_to_cluster.get(strip_ecistem_suffix(t))
        if cl is None: continue
        sc, s, e = genes[q]
        by_scaffold[sc].append((s, e, cl, bits))

with open("{output.tsv}","w") as out:
    out.write("sample\tscaffold\tstart\tend\tlength\tn_hits\tdistinct_clusters\ttop_hit_bits\tclusters_present\n")
    for sc, hits in by_scaffold.items():
        hits.sort(key=lambda x: x[0])
        cluster_buf, last_end, buckets = [], None, []
        for h in hits:
            if not cluster_buf:
                cluster_buf = [h]; last_end = h[1]; continue
            if h[0] - last_end <= MAX_GAP_BP:
                cluster_buf.append(h); last_end = max(last_end, h[1])
            else:
                buckets.append(cluster_buf); cluster_buf = [h]; last_end = h[1]
        if cluster_buf: buckets.append(cluster_buf)
        for c in buckets:
            distinct = {{x[2] for x in c}}
            if len(c) < MIN_HITS or len(distinct) < MIN_DISTINCT: continue
            cs = min(x[0] for x in c); ce = max(x[1] for x in c)
            top = max(x[3] for x in c)
            out.write("\t".join(["{params.sample}", sc, str(cs), str(ce), str(ce-cs),
                                  str(len(c)), str(len(distinct)), f"{{top:.1f}}",
                                  ",".join(sorted(distinct))]) + "\n")
PY
        """


rule ecis_collate:
    """Tier-2 deliverable: run-level eCIS candidate regions (`ecis/ecis_candidates.tsv`)."""
    input:
        expand(str(OUT/"ecis"/"{s}"/"candidates.tsv"), s=SAMPLES),
    output:
        OUT/"ecis"/"ecis_candidates.tsv",
    params:
        samples=SAMPLES,
        add_sample=False,
    conda: "../envs/python.yaml"
    script:
        "../scripts/collate_tables.py"
