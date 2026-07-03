"""Module: classify — apply v0.1.1 XGBoost models."""
rule classify_predict:
    benchmark:
        OUT/"benchmarks"/"classify_predict.tsv"
    input:
        features = OUT/"features"/"features.parquet",
    output:
        predictions = OUT/"classify"/"predictions.tsv",
    conda: "../envs/python.yaml"
    threads: 4
    script:
        "../scripts/classify_predict.py"
