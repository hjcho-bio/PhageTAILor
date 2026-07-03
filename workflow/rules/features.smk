"""Module: features — assemble per-region features for the classifier.

Inputs the six-detector union: geNomad provirus regions, tail-gene clusters, 
PHROGs-tail clusters, SecReT6 clusters, eCIStem clusters, and HMM-tail clusters.
The 79-feature schema matches Part A's training feature matrix exactly, including
the region- and genome-level HMM-derived feature blocks.
"""

rule features_assemble:
    benchmark:
        OUT/"benchmarks"/"features_assemble.tsv"
    input:
        annotations = expand(str(OUT/"annotation"/"{s}.faa"), s=SAMPLES),
        gff         = expand(str(OUT/"annotation"/"{s}.gff"), s=SAMPLES),
        secret6     = expand(str(OUT/"annotation"/"{s}.secret6.m8"), s=SAMPLES),
        ecistem     = expand(str(OUT/"annotation"/"{s}.ecistem.m8"), s=SAMPLES),
        phrogs      = expand(str(OUT/"annotation"/"{s}.phrogs.m8"), s=SAMPLES),
        tailgenes   = expand(str(OUT/"annotation"/"{s}.tailgenes.m8"), s=SAMPLES),
        hmmtail     = expand(str(OUT/"annotation"/"{s}.hmmtail.tbl"), s=SAMPLES),
        phage       = expand(str(OUT/"phage"/"{s}"/"_done"), s=SAMPLES),
        t6ss        = expand(str(OUT/"t6ss"/"{s}"/"_done"), s=SAMPLES),
        ecis        = expand(str(OUT/"ecis"/"{s}"/"_done"), s=SAMPLES),
        tailocin    = expand(str(OUT/"tailocin"/"{s}"/"_done"), s=SAMPLES),
        hmm_tail    = expand(str(OUT/"hmm_tail"/"{s}"/"_done"), s=SAMPLES),
        ani         = OUT/"taxonomy"/"nearest_train_ani.tsv",
    output:
        parquet = OUT/"features"/"features.parquet",
        regions = OUT/"features"/"regions.bed",
    conda: "../envs/python.yaml"
    threads: 4
    script:
        "../scripts/build_features.py"
