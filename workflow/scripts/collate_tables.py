"""Collate per-sample detector tables into a single run-level table.

Used by the per-detector `*_collate` rules to promote the raw per-genome
detector outputs (which live under `<detector>/{sample}/...`) into one
documented, run-level table at the output root. Empty or missing per-sample
files are skipped silently so the collation also works when a detector found
nothing in some genomes.

Snakemake wiring (see rules/*.smk):
  input        : per-sample TSVs, expanded in `SAMPLES` order
  output[0]    : collated TSV
  params.samples   : list of sample IDs, aligned 1:1 with `input`
  params.add_sample: if True, prepend a `sample` column derived from
                     params.samples (used for geNomad, whose native
                     virus_summary.tsv has no sample column); if False the
                     per-sample table is assumed to already carry `sample`.
"""
import csv
from pathlib import Path

inputs = list(snakemake.input)                       # noqa: F821  (snakemake-injected)
samples = list(snakemake.params.samples)             # noqa: F821
add_sample = bool(snakemake.params.add_sample)       # noqa: F821
out_path = Path(snakemake.output[0])                 # noqa: F821

header = None
rows = []
n_files = 0
for path, samp in zip(inputs, samples):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        continue
    with open(p) as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            h = next(reader)
        except StopIteration:
            continue
        if add_sample:
            h = ["sample"] + h
        if header is None:
            header = h
        body = list(reader)
        if not body:
            continue
        n_files += 1
        for row in body:
            if not row:
                continue
            rows.append(([samp] + row) if add_sample else row)

out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    if header is not None:
        w.writerow(header)
    w.writerows(rows)

print(f"[collate] {out_path.name}: {len(rows)} regions from {n_files}/{len(inputs)} genomes")
