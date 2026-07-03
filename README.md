# PhageTAILor v0.1.0

Predict prophages, **tailocins**, T6SS, and eCIS in bacterial genomes.

---

## Quick start

PhageTAILor ships with a beginner-friendly `phagetailor` command. You do **not**
need to know YAML, Snakemake, or MMseqs2 to run it.

### Step 1 — put your bacterial genome nucleotide FASTA files in a folder

```
my_genomes/
├── strain_A.fna
├── strain_B.fna
└── strain_C.fna
```

(`.fasta` is also accepted. Multi-contig assemblies are fine.)

### Step 2 — run

```bash
phagetailor run my_genomes/  --out my_results/  --cores 8
```

That is it. The pipeline will run the nine steps (gene calling → taxonomy →
candidate detection → classification → report) and write everything to
`my_results/`.

### Step 3 — read the results

```bash
firefox my_results/report.html       # browser-friendly summary
open    my_results/per_genome.tsv    # one row per genome (Excel-friendly)
cat     my_results/predictions.tsv   # per-region details
```

### Other input formats

If you would rather not put files in a folder, you can pass a list:

- **Tab-separated (Excel-friendly)** — a file with a `fna_path` column (and
  optionally a `sample` column to override the auto-derived sample name):

  ```
  sample    fna_path
  GCF1      /path/to/genome1.fna
  GCF2      /path/to/genome2.fna
  ```

- **Comma-separated** — same format with commas instead of tabs.

- **Plain text** — one FASTA path per line.

Run `phagetailor explain input` to see all accepted formats with examples.

### Helpful commands

```bash
phagetailor --help                   # top-level help
phagetailor example                  # a guided quick-start
phagetailor check                    # verify dependencies + databases installed
phagetailor explain annotation       # plain-language explanation of any step
phagetailor explain tailocin
phagetailor explain list             # show all explanation topics
```

Every pipeline step has a dedicated plain-language explanation accessible
via `phagetailor explain <step>`.

---

## What the pipeline does (Part B — user inference)

```
your .fna file(s)
   │
   ├─ annotation : Pyrodigal-gv → {sample}.faa + .gff
   ├─ taxonomy   : GTDB-Tk classify + skani vs training sketch
   │              → genus/family/order IDs + nearest_other_ani confidence band
   ├─ phage      : geNomad end-to-end (provirus regions)
   ├─ tailocin   : (i) MMseqs2 vs tail-gene DB → tail-gene detector
   │              (ii) MMseqs2 vs PHROGs "tail" category → PHROGs-tail detector
   │              (iii) geNomad provirus regions
   │              All three feed a candidate union; redundant detection lets us
   │              catch tailocins missing canonical phage markers.
   ├─ t6ss       : MMseqs2 vs SecReT6 → per-scaffold T6SS detector
   ├─ ecis       : MMseqs2 vs eCIStem → per-scaffold eCIS detector
   ├─ features   : 74-feature vector per candidate region (phylogenetic-placement
   │              features removed to eliminate train/deploy leakage)
   ├─ classify   : 4 LightGBM heads (multiclass + tailocin/T6SS/eCIS binary)
   └─ report     : per-genome TSV + browser-friendly HTML summary
                   (calls from binary heads @0.5 for tailocin/T6SS/eCIS)
```

Run `phagetailor explain <step>` for a beginner-friendly explanation of any step.

## Outputs

### Classifier results (produced when the `classify`/`report` modules run)

| File | What it is |
|---|---|
| `report/per_genome_predictions.tsv` | One row per input genome with max per-class scores + confidence band. **Open in Excel.** |
| `classify/predictions.tsv` | Per-region predictions: every candidate region with all class probabilities. **For genome-level filtering.** |
| `features/regions.bed` | BED file of all candidate regions for **IGV / JBrowse / UCSC** visualization. |
| `report/summary.html` | Self-contained browser-friendly summary. |
| `manifest.json` | Run provenance: tool versions, input files, model version. |

### Raw detector tables (produced by each detector module, no classifier required)

Each detector also writes one documented, run-level table collating its raw
candidate regions across all input genomes. These are first-class outputs: you
can run a single detector module and use its table directly, without invoking
the LightGBM classifier (see [Running only the steps you need](#running-only-the-steps-you-need)).

| File | Detector / module | Columns |
|---|---|---|
| `phage/prophages.tsv` | geNomad prophages (`phage`) | `sample` + geNomad `virus_summary` (seq_name, coordinates, length, n_genes, virus_score, n_hallmarks, marker_enrichment, taxonomy, …) |
| `t6ss/t6ss_candidates.tsv` | SecReT6 T6SS (`t6ss`) | sample, scaffold, start, end, length, n_hits, distinct_components, top_hit_bits, components_present |
| `ecis/ecis_candidates.tsv` | eCIStem eCIS (`ecis`) | sample, scaffold, start, end, length, n_hits, distinct_clusters, top_hit_bits, clusters_present |
| `tailocin/tailgenes_candidates.tsv` | tail-gene tailocin (`tailocin`) | sample, scaffold, start, end, length, n_hits, distinct_seqs, top_hit_bits |
| `tailocin/phrog_tail_candidates.tsv` | PHROGS-tail tailocin (`tailocin`) | sample, scaffold, start, end, length, n_hits, distinct_phrogs, top_hit_bits |
| `hmm_tail/hmm_tail_candidates.tsv` | HMM tail-cassette (`hmm_tail`) | sample, scaffold, start, end, length, n_hits, n_distinct, top_bitscore |

> The full per-genome geNomad output is still preserved under
> `phage/{sample}/genomad/` if you need every geNomad intermediate; `phage/prophages.tsv`
> is the convenient run-level summary.

## Running only the steps you need

PhageTAILor is modular: the `modules` setting selects which stages run, and
dependencies are resolved automatically. Use it to run a single detector
(e.g. prophages only) without the classifier, or to stop after candidate
detection.

| Goal | Setting |
|---|---|
| Everything (default) | `modules: "all"` |
| **Prophages only** (geNomad → `phage/prophages.tsv`) | `modules: "phage"` |
| **Tailocin candidates only** | `modules: "tailocin"` |
| **T6SS candidates only** | `modules: "t6ss"` |
| **eCIS candidates only** | `modules: "ecis"` |
| Detector tables, no classifier | `modules: "phage,tailocin,t6ss,ecis,hmm_tail"` |
| Full classification | `modules: "all"` (or `"report"`, which pulls in everything it needs) |

Set it via the `phagetailor` CLI, in `config/config.yaml`, or on the Snakemake
command line:

```bash
# Prophages only — runs geNomad and nothing else (beginner CLI)
phagetailor run samples.tsv --out results/ --cores 16 --modules phage

# Equivalent, calling Snakemake directly
snakemake -s workflow/Snakefile --cores 16 --use-conda \
    --config genome_list=my_genomes.txt output_dir=results/ modules="phage"
```

Requesting a downstream module automatically pulls in everything it depends on
(e.g. `modules: "classify"` runs all detectors → features → classifier), so you
never have to list prerequisites by hand. Conversely, every detector module
works standalone, so the only step you can't currently *omit* is an individual
detector while still running the classifier — the model was trained on all six
detectors together, so it requires the full candidate set (see the modularity
note in the project docs).

## Confidence policy

Each prediction is tagged with a confidence band reflecting how close your
query genome is to the v0.1.0 training pool:

| Band | nearest_other_ani | Interpretation |
|---|---|---|
| `high` | ≥ 95 | Test genome is close to training distribution; trust the score. |
| `medium` | 90–95 | Same genus / different species; borderline. |
| `low` | < 90 | Cross-genus extrapolation; treat scores as ranked candidates. |
| `no_match` | no aligned training genome ≥ 50% query coverage | Truly out-of-distribution. |

PhageTAILor v0.1.0 is Gammaproteobacteria-heavy (84% of training); cross-clade
predictions on, say, Cyanobacteria should be treated as low-confidence
candidates pending experimental follow-up.

## Installing

```bash
# 1. clone the repo
git clone https://github.com/hjcho-bio/PhageTAILor_final.git
cd PhageTAILor_final

# 2. create the conda environment
conda env create -f workflow/envs/python.yaml -n phagetailor
conda activate phagetailor

# 3. download the geNomad reference DB (one-time, ~5 GB)
bash scripts/setup_genomad_db.sh

# 4. download the PhageTAILor reference DBs — tail-HMM, PHROGS, and the
#    training sketch (~3.2 GB total, hosted on Zenodo)
bash scripts/setup_databases.sh

# 5. verify everything is installed
./phagetailor check
```

Reference data is distributed outside git because of its size: the geNomad DB is
fetched from geNomad, and the PhageTAILor DBs (tail-HMM library, PHROGS DB, and
training sketch) are hosted on Zenodo and fetched by `setup_databases.sh`. The
smaller detector DBs (SecReT6, tail-gene, eCIStem) ship in this repo.

The training sketch powers the confidence band (`nearest_other_ani`); you do
**not** need the 6,501 training genomes to run PhageTAILor. Their accessions are
listed in
[`resources/training_sketch/training_accessions.txt`](resources/training_sketch/training_accessions.txt)
for provenance, and `scripts/build_training_sketch.sh` can re-sketch them if you
obtain the genomes yourself.

## What if I already have annotations / geNomad output?

The pipeline can skip rules whose outputs already exist. See `config/config.yaml`
for the `pre_annotated_dir` and `pre_genomad_dir` options, or run
`phagetailor explain features` for the modular skip logic.

## Architecture (Part A — model construction)

The model was built from a training-set curation across 6,893 genomes spanning
23 phyla, a 5-detector candidate union, LightGBM training, and held-out
evaluation. The trained model files live in
[`resources/model/`](resources/model/) (`multiclass.txt`, `tailocin_binary.txt`,
`t6ss_binary.txt`, `ecis_binary.txt`) plus the encoders and feature column order;
see [`resources/model/EVAL.md`](resources/model/EVAL.md) for the full
cross-validation, algorithm-comparison, and tool-comparison results. The
training-set construction inputs, build scripts, and model-building diagnostics
are kept as supplemental data outside this repository (not needed to run the
tool).

## Data availability

The reference DBs and the training sketch are archived on Zenodo:
**[doi:10.5281/zenodo.21152308](https://doi.org/10.5281/zenodo.21152308)**.
`scripts/setup_databases.sh` downloads them automatically.

## Citing

If you use PhageTAILor, please cite both the tool and its upstream methods:

- **PhageTAILor** — this work; [Update the citation information here]
- geNomad — Camargo *et al.* 2023, *Nat. Biotechnol*, doi:10.1038/s41587-023-01953-y
- SecReT6 — Zhang *et al.* 2022, *SCIENCE CHINA Life Sciences* (and updates), doi:10.1007/s11427-022-2172-x
- eCIStem — Geller *et al.* 2021, *Nat. Commun.*, doi:10.1038/s41467-021-23777-7 (database: https://ecistem.pythonanywhere.com/)
- PHROGs — Terzian *et al.* 2021, *NAR GAB*, doi:10.1093/nargab/lqab067
- GTDB-Tk — Chaumeil *et al.* 2020, *Bioinformatics*, doi:10.1093/bioinformatics/btz848
- MMseqs2 — Steinegger & Söding 2017, *Nat. Biotechnol*, doi:10.1038/nbt.3988
- skani — Shaw & Yu 2023, *Nat. Methods*, doi:10.1038/s41592-023-02018-3
- Pyrodigal-gv — Camargo *et al.* 2023, *Nat. Biotechnol*, doi:10.1038/s41587-023-01953-y / Larralde *et al.* 2022, *JOSS*, doi:10.21105/joss.04296
- LightGBM — Ke *et al.* 2017, *NeurIPS*, doi:10.5555/3294996.3295074

## Roadmap

- Expand the held-out eCIS benchmark (currently only 2 expected-positive strains).
- Add more taxonomically diverse training data (currently 84% Gammaproteobacteria).
- Add Docker / Singularity / Apptainer containers for fully isolated installs.
- Ship a tiny example genome so `phagetailor run examples/ --out test/` works
  out-of-the-box without users having to download anything.

## Advanced usage (Snakemake direct)

For experienced users who want fine-grained control, the underlying Snakemake
workflow lives in [`workflow/`](workflow/) and can be invoked directly:

```bash
snakemake -s workflow/Snakefile --cores 32 \
    --config genome_list=my_genomes.txt output_dir=results/ \
    --use-conda
```

A legacy bash wrapper is also preserved at `phagetailor.legacy.sh` for users
of the v0.1.1 interface.

## License

PhageTAILor is released under the MIT License — see [`LICENSE`](LICENSE).
