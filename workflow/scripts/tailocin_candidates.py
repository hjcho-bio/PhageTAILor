"""Tailocin candidate detection: stub.

For each genome:
  1. Read geNomad-detected provirus regions (from {OUT}/phage/{s}/genomad/...).
  2. For each region, score "tailocin-likeness" using:
     - hits to tail/baseplate/sheath PHROGs (TODO: bundle PHROG HMM/MMseqs lookup)
     - hits to SecReT6 (any T6SS hit suggests NOT tailocin)
     - region length (10–25 kb typical for tailocin)
     - rule-based exclusion: capsid/portal genes => not a tailocin
  3. Write candidates.bed: scaffold, start, end, name, provisional_score, strand
"""
from pathlib import Path
import sys

# ── Snakemake-injected ──────────────────────────────────────────────────────
sample = snakemake.params.sample
outdir = Path(snakemake.params.outdir)
secret6_m8 = Path(snakemake.input.secret6_hits)
gff = Path(snakemake.input.gff)
out_bed = Path(snakemake.output.bed)

outdir.mkdir(parents=True, exist_ok=True)

# ── TODO Day 2-3: implement scoring ────────────────────────────────────────
# Stub: emit an empty BED so downstream rules don't break.
out_bed.write_text("")
print(f"[tailocin_candidates] sample={sample}: stub – write features in features.smk", file=sys.stderr)
