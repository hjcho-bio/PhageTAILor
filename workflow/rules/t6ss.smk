"""Module: t6ss — cluster SecReT6 MMseqs hits per scaffold into T6SS candidate regions.

This replaces the (deferred) MacSyFinder dependency for v0.1 and matches exactly
how the model was trained (same clustering parameters).
"""

rule t6ss_cluster_secret6:
    benchmark:
        OUT/"benchmarks"/"t6ss_cluster_secret6.{s}.tsv"
    input:
        m8  = OUT/"annotation"/"{s}.secret6.m8",
        gff = OUT/"annotation"/"{s}.gff",
    output:
        flag = touch(OUT/"t6ss"/"{s}"/"_done"),
        tsv  = OUT/"t6ss"/"{s}"/"candidates.tsv",
    conda: "../envs/python.yaml"
    params:
        meta = config["secret6_metadata"],
        sample = lambda wc: wc.s,
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        python3 - <<PY
import csv, re
from collections import defaultdict
from pathlib import Path

CORE = ["TssA","TssB","TssC","TssD","TssE","TssF","TssG","TssH","TssI","TssJ","TssK","TssL","TssM"]
MIN_BITS, MIN_PIDENT = 50, 30
MAX_GAP_BP = 15000
MIN_DISTINCT, MIN_HITS = 3, 4

# component map
target_to_comp = {{}}
with open("{params.meta}") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd: target_to_comp[r["protein_id"]] = r["component"]

# gene coords from GFF
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

# read m8 hits, top per query
seen = set()
by_scaffold = defaultdict(list)
with open("{input.m8}") as fh:
    for line in fh:
        f = line.rstrip().split("\t")
        if len(f) < 6: continue
        q, t = f[0], f[1]
        try: pid, alnlen, evalue, bits = float(f[2]), int(f[3]), float(f[4]), float(f[5])
        except ValueError: continue
        if bits < MIN_BITS or pid < MIN_PIDENT: continue
        if q in seen: continue
        seen.add(q)
        if q not in genes: continue
        sc, s, e = genes[q]
        comp = target_to_comp.get(t, "?")
        by_scaffold[sc].append((s, e, comp, bits))

with open("{output.tsv}","w") as out:
    out.write("sample\tscaffold\tstart\tend\tlength\tn_hits\tdistinct_components\ttop_hit_bits\tcomponents_present\n")
    for sc, hits in by_scaffold.items():
        hits.sort(key=lambda x: x[0])
        cluster, last_end, clusters = [], None, []
        for h in hits:
            if not cluster:
                cluster = [h]; last_end = h[1]; continue
            if h[0] - last_end <= MAX_GAP_BP:
                cluster.append(h); last_end = max(last_end, h[1])
            else:
                clusters.append(cluster); cluster = [h]; last_end = h[1]
        if cluster: clusters.append(cluster)
        for c in clusters:
            distinct = {{x[2] for x in c}} - {{"?"}}
            if len(c) < MIN_HITS: continue
            if len(distinct) < MIN_DISTINCT: continue
            cs = min(x[0] for x in c); ce = max(x[1] for x in c)
            top = max(x[3] for x in c)
            out.write("\t".join(["{params.sample}", sc, str(cs), str(ce), str(ce-cs),
                                  str(len(c)), str(len(distinct)), f"{{top:.1f}}",
                                  ",".join(sorted(distinct))]) + "\n")
PY
        """


rule t6ss_collate:
    """Tier-2 deliverable: run-level T6SS candidate regions (`t6ss/t6ss_candidates.tsv`)."""
    input:
        expand(str(OUT/"t6ss"/"{s}"/"candidates.tsv"), s=SAMPLES),
    output:
        OUT/"t6ss"/"t6ss_candidates.tsv",
    params:
        samples=SAMPLES,
        add_sample=False,
    conda: "../envs/python.yaml"
    script:
        "../scripts/collate_tables.py"
