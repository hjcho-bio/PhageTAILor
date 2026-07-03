"""Module: phage — geNomad end-to-end. Produces virus_summary.tsv (the candidate source).

Inputs are symlinked to {sample}.fna so geNomad's output directory is named consistently.
"""

rule genomad_run:
    benchmark:
        OUT/"benchmarks"/"genomad_run.{s}.tsv"
    input:
        fna = lambda wc: GENOMES[wc.s],
    output:
        flag    = touch(OUT/"phage"/"{s}"/"_done"),
        summary = OUT/"phage"/"{s}"/"genomad"/"{s}_summary"/"{s}_virus_summary.tsv",
    threads: config["threads"]["genomad"]
    conda: "../envs/phage.yaml"
    params:
        db      = config["genomad_db"],
        relaxed = "--relaxed" if config.get("genomad_relaxed", True) else "",
        outdir  = lambda wc: str(OUT/"phage"/wc.s/"genomad"),
        sample  = lambda wc: wc.s,
    shell:
        r"""
        mkdir -p {params.outdir}
        # Symlink the input under a clean basename so geNomad's output paths are stable
        clean_fna={params.outdir}/{params.sample}.fna
        ln -sf $(readlink -f {input.fna}) ${{clean_fna}}
        genomad end-to-end {params.relaxed} --threads {threads} \
                --cleanup \
                ${{clean_fna}} {params.outdir} {params.db}
        """


rule phage_collate:
    """Tier-1 deliverable: promote geNomad's per-genome virus_summary.tsv into one
    documented run-level prophage table (`phage/prophages.tsv`). Available whenever
    the `phage` module runs, including prophage-only runs (`modules: annotation,phage`)
    that never touch the classifier."""
    input:
        expand(str(OUT/"phage"/"{s}"/"genomad"/"{s}_summary"/"{s}_virus_summary.tsv"), s=SAMPLES),
    output:
        OUT/"phage"/"prophages.tsv",
    params:
        samples=SAMPLES,
        add_sample=True,
    conda: "../envs/python.yaml"
    script:
        "../scripts/collate_tables.py"
