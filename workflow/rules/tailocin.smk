"""Module: tailocin — cluster tail-gene MMseqs hits and PHROGs-tail hits per scaffold.

Two independent cluster-based detectors feed tailocin candidates:
  1. tail-gene clusterer — MMseqs hits vs the bundled training-tailocin protein DB
     (`resources/tailgenes_db/`). HMM-free; mirrors the training-side detector.
  2. PHROGs-tail clusterer — MMseqs hits vs PHROGs, restricted to category="tail".

Both are independent of geNomad and recover tailocins that lack canonical phage
hallmarks (the dominant failure mode of geNomad-only tailocin detection).
The geNomad-detected provirus regions are also fed into the candidate union
via the `phage` rule, so the tailocin candidate set is the union of THREE
detectors at inference time.
"""

rule tailocin_cluster_tailgenes:
    benchmark:
        OUT/"benchmarks"/"tailocin_cluster_tailgenes.{s}.tsv"
    input:
        m8  = OUT/"annotation"/"{s}.tailgenes.m8",
        gff = OUT/"annotation"/"{s}.gff",
    output:
        tsv = OUT/"tailocin"/"{s}"/"tailgenes_candidates.tsv",
    conda: "../envs/python.yaml"
    params:
        sample = lambda wc: wc.s,
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        python3 - <<PY
import re
from collections import defaultdict

MIN_BITS, MIN_PIDENT = 50, 30
MAX_GAP_BP = 15000
MIN_DISTINCT, MIN_HITS = 3, 4

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
        seqid = t.split("|", 1)[0]
        sc, s, e = genes[q]
        by_scaffold[sc].append((s, e, seqid, bits))

with open("{output.tsv}","w") as out:
    out.write("sample\tscaffold\tstart\tend\tlength\tn_hits\tdistinct_seqs\ttop_hit_bits\n")
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
                                  str(len(c)), str(len(distinct)), f"{{top:.1f}}"]) + "\n")
PY
        """

rule tailocin_cluster_phrogs:
    benchmark:
        OUT/"benchmarks"/"tailocin_cluster_phrogs.{s}.tsv"
    input:
        m8  = OUT/"annotation"/"{s}.phrogs.m8",
        gff = OUT/"annotation"/"{s}.gff",
    output:
        tsv = OUT/"tailocin"/"{s}"/"phrog_tail_candidates.tsv",
    conda: "../envs/python.yaml"
    params:
        annot = config["phrogs_annotation"],
        sample = lambda wc: wc.s,
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        python3 - <<PY
import csv, re
from collections import defaultdict

MIN_BITS = 50
MAX_GAP_BP = 15000
MIN_DISTINCT_TAIL, MIN_HITS_TAIL = 2, 3   # looser thresholds: tail-gene cassettes are short

# PHROG ID → category; we only keep "tail"
phrog_cat = {{}}
with open("{params.annot}") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        try: phrog_cat[int(r["phrog"])] = r.get("category","")
        except: pass

def parse_phrog(t):
    if t.startswith("phrog_"):
        try: return int(t[len("phrog_"):])
        except: return None
    return None

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
        try: bits = float(f[5])
        except ValueError: continue
        if bits < MIN_BITS: continue
        if q in seen: continue
        seen.add(q)
        if q not in genes: continue
        ph = parse_phrog(t)
        if ph is None: continue
        if phrog_cat.get(ph,"") != "tail": continue
        sc, s, e = genes[q]
        by_scaffold[sc].append((s, e, ph, bits))

with open("{output.tsv}","w") as out:
    out.write("sample\tscaffold\tstart\tend\tlength\tn_hits\tdistinct_phrogs\ttop_hit_bits\n")
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
            if len(c) < MIN_HITS_TAIL or len(distinct) < MIN_DISTINCT_TAIL: continue
            cs = min(x[0] for x in c); ce = max(x[1] for x in c)
            top = max(x[3] for x in c)
            out.write("\t".join(["{params.sample}", sc, str(cs), str(ce), str(ce-cs),
                                  str(len(c)), str(len(distinct)), f"{{top:.1f}}"]) + "\n")
PY
        """

rule tailocin_marker:
    benchmark:
        OUT/"benchmarks"/"tailocin_marker.{s}.tsv"
    input:
        OUT/"phage"/"{s}"/"_done",
        OUT/"tailocin"/"{s}"/"tailgenes_candidates.tsv",
        OUT/"tailocin"/"{s}"/"phrog_tail_candidates.tsv",
    output:
        touch(OUT/"tailocin"/"{s}"/"_done"),
    shell:
        r"""mkdir -p $(dirname {output[0]}); :"""


rule tailocin_collate_tailgenes:
    """Tier-2 deliverable: run-level tail-gene candidate regions
    (`tailocin/tailgenes_candidates.tsv`)."""
    input:
        expand(str(OUT/"tailocin"/"{s}"/"tailgenes_candidates.tsv"), s=SAMPLES),
    output:
        OUT/"tailocin"/"tailgenes_candidates.tsv",
    params:
        samples=SAMPLES,
        add_sample=False,
    conda: "../envs/python.yaml"
    script:
        "../scripts/collate_tables.py"


rule tailocin_collate_phrogs:
    """Tier-2 deliverable: run-level PHROGS-tail candidate regions
    (`tailocin/phrog_tail_candidates.tsv`)."""
    input:
        expand(str(OUT/"tailocin"/"{s}"/"phrog_tail_candidates.tsv"), s=SAMPLES),
    output:
        OUT/"tailocin"/"phrog_tail_candidates.tsv",
    params:
        samples=SAMPLES,
        add_sample=False,
    conda: "../envs/python.yaml"
    script:
        "../scripts/collate_tables.py"
