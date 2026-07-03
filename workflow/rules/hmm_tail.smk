"""Module: hmm_tail — cluster HMM tail-cassette hits per scaffold.

This is the HMM (HMMER profile) counterpart to the tail-gene MMseqs clusterer.
Hits to the 312-profile cassette-only detector subset (baseplate, sheath, tube,
fibre, holin, lysin) are merged into candidate regions using the same
per-scaffold gap-merging rule as the PHROGs-tail clusterer:
  - top-bitscore HMM hit per query gene (after E-value filtering)
  - sort by start coordinate; merge consecutive hits whose gene-to-gene gap is ≤ 15 kb
  - admit a merged window with ≥ 3 hits spanning ≥ 2 distinct HMM profiles

The HMM clusterer is the *sequence-divergence-tolerant* counterpart to the
sequence-based detectors and is the primary mechanism for picking up cross-clade
tail loci where MMseqs2 sequence search falls below threshold.
"""

rule hmm_tail_cluster:
    benchmark:
        OUT/"benchmarks"/"hmm_tail_cluster.{s}.tsv"
    input:
        tbl = OUT/"annotation"/"{s}.hmmtail.tbl",
        gff = OUT/"annotation"/"{s}.gff",
    output:
        tsv  = OUT/"hmm_tail"/"{s}"/"hmm_tail_candidates.tsv",
        done = touch(OUT/"hmm_tail"/"{s}"/"_done"),
    conda: "../envs/python.yaml"
    params:
        sample      = lambda wc: wc.s,
        hmm_meta    = config["hmm_tail_metadata"],   # detector-subset metadata (312 HMMs)
        max_gap_bp  = 15000,
        min_hits    = 3,
        min_distinct= 2,
    script:
        "../scripts/hmm_tail_cluster.py"


rule hmm_tail_collate:
    """Tier-2 deliverable: run-level HMM tail-cassette candidate regions
    (`hmm_tail/hmm_tail_candidates.tsv`)."""
    input:
        expand(str(OUT/"hmm_tail"/"{s}"/"hmm_tail_candidates.tsv"), s=SAMPLES),
    output:
        OUT/"hmm_tail"/"hmm_tail_candidates.tsv",
    params:
        samples=SAMPLES,
        add_sample=False,
    conda: "../envs/python.yaml"
    script:
        "../scripts/collate_tables.py"
