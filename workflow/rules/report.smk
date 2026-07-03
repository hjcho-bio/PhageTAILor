"""Module: report — per-genome TSV + summary HTML."""
rule report_summary:
    benchmark:
        OUT/"benchmarks"/"report_summary.tsv"
    input:
        predictions = OUT/"classify"/"predictions.tsv",
        regions     = OUT/"features"/"regions.bed",
    output:
        tsv  = OUT/"report"/"per_genome_predictions.tsv",
        html = OUT/"report"/"summary.html",
    conda: "../envs/python.yaml"
    script:
        "../scripts/build_report.py"
