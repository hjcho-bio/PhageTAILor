"""Per-region feature assembly for inference (Part B).

Mirrors the training-side feature builder so the inference feature matrix
matches the schema the model was trained on (79 features).

Inputs (via snakemake.input):
  annotations  : list of {sample}.faa paths
  gff          : list of {sample}.gff paths
  secret6      : list of {sample}.secret6.m8 paths
  ecistem      : list of {sample}.ecistem.m8 paths
  phrogs       : list of {sample}.phrogs.m8 paths
  tailgenes    : list of {sample}.tailgenes.m8 paths
  hmmtail      : list of {sample}.hmmtail.tbl paths      (HMMER tblout)
  ani          : taxonomy/nearest_train_ani.tsv

Detector union (six): geNomad provirus + tail-gene cluster + PHROGs-tail cluster + 
SecReT6 cluster + eCIStem cluster + HMM-tail cluster.
"""
import csv, json, os, re, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd

OUT_PARQUET = Path(snakemake.output.parquet)
OUT_REGIONS = Path(snakemake.output.regions)
OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)

CFG = snakemake.config
RUN_DIR = Path(CFG["output_dir"])
SAMPLES = sorted({Path(p).stem for p in snakemake.input.annotations})

CORE_TSS = ["TssA","TssB","TssC","TssD","TssE","TssF","TssG","TssH","TssI","TssJ","TssK","TssL","TssM"]
ACC_TSS  = ["PAAR","FHA","TagF","VgrG","Hcp","Tae4","TssN","TssO"]
ALL_COMP = CORE_TSS + ACC_TSS
PHROG_CATS = ["tail","head and packaging","DNA, RNA and nucleotide metabolism",
              "lysis","integration and excision","transcription regulation",
              "moron, auxiliary metabolic gene and host takeover","connector","other","unknown function"]
MIN_BITS, MIN_PIDENT = 50.0, 30.0
MAX_GAP_BP = 15000
MIN_DISTINCT, MIN_HITS = 3, 4

def canon(s): return re.sub(r'\.\d+$', '', s)

# ─── Load reference metadata ──────────────────────────────────────────────
sec6_comp = {}
with open(CFG["secret6_metadata"]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sec6_comp[r["protein_id"]] = r["component"]

ecis_clust = {}
with open(CFG["ecistem_metadata"]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ecis_clust[r["seqid"]] = r["cluster_id"]

phrog_cat = {}
with open(CFG["phrogs_annotation"]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        try: phrog_cat[int(r["phrog"])] = r.get("category","unknown function")
        except: pass

# HMM tail-cassette metadata (HMM name -> module category) for the full
# 597-profile library (full module set, including capsid/portal as
# discriminative features). Detector-subset metadata is used inside the
# HMM-tail clusterer rule, not here.
hmm_cat = {}
HMM_DET_CATS  = {"baseplate","sheath","tube","fibre","holin","lysin","tailocin_spec","tail_other"}
HMM_REGION_MODS = ["baseplate","sheath","tube","fibre","lysin"]
HMM_GENOME_MODS = ["baseplate","sheath","tube","fibre","holin","lysin","capsid","portal"]
HMM_FULL_META = CFG.get("hmm_tail_full_metadata") or "resources/hmm_tail_db/hmm_full_tail.tsv"
if Path(HMM_FULL_META).exists():
    with open(HMM_FULL_META) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            hmm_cat[r["hmm"]] = r.get("category","")

# ─── Load taxonomy + ANI ──────────────────────────────────────────────────
sample_tax = {}
tax_path = RUN_DIR/"taxonomy"/"gtdbtk.bac120.summary.tsv"
if tax_path.exists() and tax_path.stat().st_size > 0:
    with open(tax_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ug = r.get("user_genome",""); cls = r.get("classification","")
            parts = {}
            for p in cls.split(";"):
                if "__" in p:
                    k, v = p.split("__",1); parts[k.strip()] = v.strip()
            sample_tax[ug] = (parts.get("g","unknown"), parts.get("f","unknown"), parts.get("o","unknown"))

ani_path = RUN_DIR/"taxonomy"/"nearest_train_ani.tsv"
sample_ani = {}; sample_n_close = defaultdict(int)
if ani_path.exists():
    with open(ani_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sample_ani[r["sample"]] = float(r.get("nearest_other_ani") or 0.0)
            sample_n_close[r["sample"]] = int(r.get("n_close_train") or 0)

with open(Path(CFG["model_dir"])/"encoders.json") as fh:
    ENC = json.load(fh)

# ─── Helpers ──────────────────────────────────────────────────────────────
def parse_phrog_target(t):
    if t.startswith("phrog_"):
        try: return int(t[len("phrog_"):])
        except: return None
    return None
def strip_ecistem_suffix(t): return t.split("|", 1)[0]

def gene_index_from_gff(sample):
    """Returns {gene_id: (scaffold, start, end, strand)}. strand: 1=+, -1=-, 0=unknown."""
    gff = RUN_DIR/"annotation"/f"{sample}.gff"
    out = {}
    if not gff.exists(): return out
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS": continue
            attrs = dict(p.split("=",1) for p in f[8].split(";") if "=" in p)
            gid = attrs.get("ID","").rstrip(";")
            if not gid: continue
            try: s, e = int(f[3]), int(f[4])
            except ValueError: continue
            st = f[6].strip()
            strand = 1 if st == "+" else (-1 if st == "-" else 0)
            out[gid] = (canon(f[0]), min(s,e), max(s,e), strand)
    return out

def load_top_hits(m8_path, target_xform, min_pident=MIN_PIDENT):
    out = {}
    if not Path(m8_path).exists(): return out
    with open(m8_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6: continue
            q, t = f[0], f[1]
            try: pid, bits = float(f[2]), float(f[5])
            except ValueError: continue
            if bits < MIN_BITS or pid < min_pident: continue
            v = target_xform(t, bits)
            if v is None: continue
            prev = out.get(q)
            if prev is None or bits > prev[-1]: out[q] = v
    return out

def load_hmm_top(tbl_path, min_bits=30.0):
    """Parse hmmsearch --tblout: per-gene top HMM hit (any of the 597 profiles).
    Returns {gene_id: (hmm_name, score, category)}."""
    out = {}
    if not Path(tbl_path).exists(): return out
    with open(tbl_path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.split()
            if len(f) < 6: continue
            gene = f[0]; hmm = f[2]
            try: score = float(f[5])
            except (IndexError, ValueError): continue
            if score < min_bits: continue
            cat = hmm_cat.get(hmm, "")
            prev = out.get(gene)
            if prev is None or score > prev[1]:
                out[gene] = (hmm, score, cat)
    return out

def load_hmm_tail_regions(sample):
    """Load HMM-tail clusterer regions (the 6th detector)."""
    p = RUN_DIR/"hmm_tail"/sample/"hmm_tail_candidates.tsv"
    out = []
    if not p.exists(): return out
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                out.append({"scaffold": canon(r["scaffold"]),
                            "start": int(r["start"]), "end": int(r["end"]),
                            "length": int(r["length"]), "n_hits": int(r["n_hits"]),
                            "n_distinct": int(r["n_distinct"]),
                            "top_bits": float(r["top_bitscore"]),
                            "src": "hmm_tail",
                            "virus_score": 0.0, "n_hallmarks": 0,
                            "marker_enrichment": 0.0, "is_caudoviricetes": 0,
                            "n_genes": int(r["n_hits"])})
            except (ValueError, KeyError):
                continue
    return out

def cluster_hits(scaffold_genes, gene_to_value, min_distinct, min_hits):
    clusters = []
    for sc, glist in scaffold_genes.items():
        ann = [(s, e, gid, gene_to_value[gid]) for s, e, gid in glist if gid in gene_to_value]
        if not ann: continue
        ann.sort()
        cur, last_end = [], None
        buckets = []
        for s, e, gid, v in ann:
            if not cur:
                cur = [(s, e, gid, v)]; last_end = e; continue
            if s - last_end <= MAX_GAP_BP:
                cur.append((s, e, gid, v)); last_end = max(last_end, e)
            else:
                buckets.append(cur); cur = [(s, e, gid, v)]; last_end = e
        if cur: buckets.append(cur)
        for c in buckets:
            distinct = {x[3][0] if isinstance(x[3], tuple) else x[3] for x in c}
            if len(distinct) < min_distinct: continue
            if len(c) < min_hits: continue
            cs = min(x[0] for x in c); ce = max(x[1] for x in c)
            top_b = max(x[3][-1] if isinstance(x[3], tuple) else 0 for x in c)
            clusters.append({"scaffold": sc, "start": cs, "end": ce,
                              "length": ce - cs, "n_hits": len(c),
                              "n_distinct": len(distinct), "top_bits": top_b,
                              "components": list(distinct)})
    return clusters

def load_genomad_candidates(sample):
    p = RUN_DIR/"phage"/sample/"genomad"/f"{sample}_summary"/f"{sample}_virus_summary.tsv"
    out = []
    if not p.exists(): return out
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq = r["seq_name"]; coord = r.get("coordinates","NA")
            if "|provirus_" in seq:
                sc, prov = seq.split("|provirus_")
                try: s, e = map(int, prov.split("_"))
                except: continue
                sc = canon(sc)
            else:
                sc = canon(seq)
                if coord and coord != "NA":
                    try: s, e = map(int, coord.split("-"))
                    except: continue
                else:
                    try: s, e = 1, int(r["length"])
                    except: continue
            out.append({"scaffold": sc, "start": s, "end": e,
                         "length": int(r.get("length", e-s)),
                         "n_genes": int(r.get("n_genes", 0)),
                         "virus_score": float(r.get("virus_score", 0) or 0),
                         "n_hallmarks": int(r.get("n_hallmarks", 0) or 0),
                         "marker_enrichment": float(r.get("marker_enrichment", 0) or 0),
                         "is_caudoviricetes": int("Caudoviricetes" in r.get("taxonomy","")),
                         "src": "genomad"})
    return out

# ─── Build per-genome rows ────────────────────────────────────────────────
all_rows = []; bed_rows = []
for s in SAMPLES:
    genes = gene_index_from_gff(s)
    by_sc = defaultdict(list)
    for gid, (sc, gs, ge, _st) in genes.items(): by_sc[sc].append((gs, ge, gid))
    for sc in by_sc: by_sc[sc].sort()

    sec6_top  = load_top_hits(RUN_DIR/"annotation"/f"{s}.secret6.m8",
                                lambda t, bits: ((sec6_comp[t], bits) if t in sec6_comp else None))
    ecis_top  = load_top_hits(RUN_DIR/"annotation"/f"{s}.ecistem.m8",
                                lambda t, bits: ((ecis_clust[strip_ecistem_suffix(t)], bits) if strip_ecistem_suffix(t) in ecis_clust else None))
    phrog_top = load_top_hits(RUN_DIR/"annotation"/f"{s}.phrogs.m8",
                                lambda t, bits: (((phrog_cat.get(parse_phrog_target(t),"unknown function"), bits, parse_phrog_target(t)) if parse_phrog_target(t) is not None else None)),
                                min_pident=0)
    tail_top  = load_top_hits(RUN_DIR/"annotation"/f"{s}.tailgenes.m8",
                                lambda t, bits: (t.split("|",1)[0], bits))
    hmm_top   = load_hmm_top(RUN_DIR/"annotation"/f"{s}.hmmtail.tbl")

    # Genome-level HMM feature block (constant across regions of this sample)
    cat_counts_g = defaultdict(int)
    distinct_det = set()
    for gid, (hmm, sc_, cat) in hmm_top.items():
        cat_counts_g[cat] += 1
        if cat in HMM_DET_CATS: distinct_det.add(hmm)
    g_hmm = {"g_hmm_distinct": len(distinct_det)}
    for m in HMM_GENOME_MODS:
        g_hmm[f"g_hmm_{m}"] = int(cat_counts_g.get(m, 0))

    g_regions = load_genomad_candidates(s)
    for c in cluster_hits(by_sc, sec6_top, MIN_DISTINCT, MIN_HITS):
        comps = set(c["components"])
        c.update({"src":"secret6","virus_score":0.0,"n_hallmarks":0,
                  "marker_enrichment":0.0,"is_caudoviricetes":0,"n_genes":c["n_hits"]})
        for comp in ALL_COMP: c[f"has_{comp}"] = int(comp in comps)
        c["core_completeness"] = sum(1 for cc in CORE_TSS if cc in comps)/len(CORE_TSS)
        g_regions.append(c)
    for c in cluster_hits(by_sc, ecis_top, MIN_DISTINCT, MIN_HITS):
        c.update({"src":"ecistem","virus_score":0.0,"n_hallmarks":0,
                  "marker_enrichment":0.0,"is_caudoviricetes":0,"n_genes":c["n_hits"]})
        g_regions.append(c)
    tail_only_phrog = {gid: (v[2], v[1]) for gid, v in phrog_top.items() if v[0] == "tail"}
    for c in cluster_hits(by_sc, tail_only_phrog, 2, 3):
        c.update({"src":"phrog_tail","virus_score":0.0,"n_hallmarks":0,
                  "marker_enrichment":0.0,"is_caudoviricetes":0,"n_genes":c["n_hits"]})
        g_regions.append(c)
    for c in cluster_hits(by_sc, tail_top, MIN_DISTINCT, MIN_HITS):
        c.update({"src":"tailgenes","virus_score":0.0,"n_hallmarks":0,
                  "marker_enrichment":0.0,"is_caudoviricetes":0,"n_genes":c["n_hits"]})
        g_regions.append(c)
    # 6th detector: HMM-tail clusterer regions (already pre-computed by rules/hmm_tail.smk)
    g_regions.extend(load_hmm_tail_regions(s))

    g, fa, o = sample_tax.get(s, ("unknown","unknown","unknown"))
    genus_id  = ENC["genus"].get(g, 0)
    family_id = ENC["family"].get(fa, 0)
    order_id  = ENC["order"].get(o, 0)
    nani = sample_ani.get(s, 0.0); n_close = sample_n_close.get(s, 0)

    for reg in g_regions:
        if "has_TssA" not in reg:
            for c in ALL_COMP: reg[f"has_{c}"] = 0
            reg["core_completeness"] = 0.0
        sc = reg["scaffold"]
        try: s0, e0 = int(reg["start"]), int(reg["end"])
        except: s0, e0 = 0, 0
        in_region = [gid for (gs, ge, gid) in by_sc.get(sc, []) if not (ge < s0 or gs > e0)]
        sec6_in = [sec6_top[gid] for gid in in_region if gid in sec6_top]
        ecis_in = [ecis_top[gid] for gid in in_region if gid in ecis_top]
        phrog_in = [phrog_top[gid] for gid in in_region if gid in phrog_top]
        tail_in = [tail_top[gid] for gid in in_region if gid in tail_top]
        sec6_distinct = {c for c, b in sec6_in}

        cat_counts = {f"n_phrog_{c.replace(', ','_').replace(' ','_')}": 0 for c in PHROG_CATS}
        for gid in in_region:
            v = phrog_top.get(gid)
            if v: cat_counts[f"n_phrog_{v[0].replace(', ','_').replace(' ','_')}"] += 1
        reg.update(cat_counts)

        reg["n_secret6_hits"] = len(sec6_in)
        reg["n_secret6_distinct"] = len(sec6_distinct)
        reg["secret6_top_bits"] = max((b for c,b in sec6_in), default=0.0)
        reg["n_ecistem_hits"] = len(ecis_in)
        reg["n_ecistem_distinct_clusters"] = len({c for c,b in ecis_in})
        reg["ecistem_top_bits"] = max((b for c,b in ecis_in), default=0.0)
        reg["phrog_top_bits"] = max((b for c,b,p in phrog_in), default=0.0)
        reg["n_tailgenes_hits"] = len(tail_in)
        reg["n_tailgenes_distinct_seqs"] = len({sid for sid,b in tail_in})
        reg["tailgenes_top_bits"] = max((b for sid,b in tail_in), default=0.0)
        reg["core_completeness_T6SS"] = sum(1 for c in CORE_TSS if c in sec6_distinct)/len(CORE_TSS)
        for c in ALL_COMP:
            if f"has_{c}" not in reg: reg[f"has_{c}"] = int(c in sec6_distinct)

        # HMM region-level features
        hmm_in = [hmm_top[gid] for gid in in_region if gid in hmm_top]
        hmm_det = [(h, sc_, cat) for (h, sc_, cat) in hmm_in if cat in HMM_DET_CATS]
        reg["hmm_n_hits"]     = len(hmm_det)
        reg["hmm_n_distinct"] = len({h for h, _, _ in hmm_det})
        reg["hmm_top_bits"]   = max((sc_ for _, sc_, _ in hmm_det), default=0.0)
        present_mods = {cat for _, _, cat in hmm_in}
        for m in HMM_REGION_MODS:
            reg[f"hmm_has_{m}"] = int(m in present_mods)
        # Genome-level HMM features (broadcast same value to every region of this sample)
        reg.update(g_hmm)

        # Gene-architecture features (n_genes_per_kb, strand_consistency,
        # tail_gene_span_frac, tail_gene_compactness, lysis_near_end).
        # Mirrors the training-side compute_arch() recipe exactly.
        length_bp = max(1, e0 - s0)
        n_in = len(in_region)
        if n_in == 0:
            reg["n_genes_per_kb"] = 0.0
            reg["strand_consistency"] = 0.5
            reg["tail_gene_span_frac"] = 0.0
            reg["tail_gene_compactness"] = 0.0
            reg["lysis_near_end"] = 0
        else:
            reg["n_genes_per_kb"] = n_in / (length_bp / 1000.0)
            plus = sum(1 for gid in in_region if gid in genes and genes[gid][3] == 1)
            minus = sum(1 for gid in in_region if gid in genes and genes[gid][3] == -1)
            reg["strand_consistency"] = max(plus, minus) / n_in
            tail_mids = [(genes[gid][1]+genes[gid][2])//2
                          for gid in in_region if gid in tail_top and gid in genes]
            if len(tail_mids) >= 2:
                sm = sorted(tail_mids)
                reg["tail_gene_span_frac"] = (sm[-1]-sm[0]) / length_bp
                gaps = [sm[i+1]-sm[i] for i in range(len(sm)-1)]
                mean_g = sum(gaps)/len(gaps)
                std_g = (sum((gp-mean_g)**2 for gp in gaps)/len(gaps))**0.5
                reg["tail_gene_compactness"] = 1.0 - std_g/length_bp
            else:
                reg["tail_gene_span_frac"] = 0.0
                reg["tail_gene_compactness"] = 0.0
            thr = s0 + 0.8*length_bp
            lne = 0
            for gid in in_region:
                v = phrog_top.get(gid)
                if v and v[0] == "lysis" and gid in genes:
                    if (genes[gid][1]+genes[gid][2])//2 >= thr:
                        lne = 1; break
            reg["lysis_near_end"] = lne

        reg["sample"] = s
        reg["genus_id"] = genus_id; reg["family_id"] = family_id; reg["order_id"] = order_id
        reg["nearest_other_ani"] = nani; reg["n_close_train"] = n_close
        reg["candidate_id"] = f"{s}:{reg['scaffold']}:{reg['start']}-{reg['end']}:{reg['src']}"
        all_rows.append(reg)
        bed_rows.append((reg["scaffold"], reg["start"], reg["end"], f"{s}|{reg['src']}"))

df = pd.DataFrame(all_rows).fillna(0)
# Coerce object columns to str so pyarrow can serialize them (fillna(0) can leave
# mixed str/int object columns, which made to_parquet fail -> a misleading 0-byte file).
for _c in df.columns:
    if df[_c].dtype == object:
        df[_c] = df[_c].astype(str)
# Always write the tsv.gz (canonical fallback that classify reads), then try parquet.
df.to_csv(str(OUT_PARQUET).replace(".parquet", ".tsv.gz"), sep="\t", index=False, compression="gzip")
try:
    df.to_parquet(OUT_PARQUET, index=False)
except Exception as _e:
    import sys as _sys
    print(f"[build_features] parquet write failed ({_e}); using tsv.gz fallback", file=_sys.stderr)
    OUT_PARQUET.write_bytes(b"")
with open(OUT_REGIONS, "w") as fh:
    for sc, st, en, name in bed_rows: fh.write(f"{sc}\t{st}\t{en}\t{name}\n")
print(f"[features] wrote {len(df)} regions across {len(SAMPLES)} samples", file=sys.stderr)
print(f"[features] sources: {df['src'].value_counts().to_dict() if len(df) else '{}'}", file=sys.stderr)
