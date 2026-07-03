"""Module: taxonomy — fast ANI to training sketch (default).

Optionally GTDB-Tk if config['use_gtdbtk'] is true (slow, days-scale).
"""

rule nearest_train_ani:
    benchmark:
        OUT/"benchmarks"/"nearest_train_ani.tsv"
    """Compute nearest-ANI to the staged training-genome sketch (resources/training_sketch).

    Produces results/taxonomy/nearest_train_ani.tsv:
       sample, nearest_other_ani, n_close_train, top_match
    """
    input:
        fnas = list(GENOMES.values()),
    output:
        tsv = OUT/"taxonomy"/"nearest_train_ani.tsv",
    threads: config["threads"]["gtdbtk"]
    conda: "../envs/gtdbtk.yaml"
    params:
        sketch = "resources/training_sketch/sketch_db",
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        tmp=$(mktemp -d); trap "rm -rf $tmp" EXIT
        if [ ! -e {params.sketch} ]; then
            echo "ERROR: training sketch not found at {params.sketch}." >&2
            echo "Download it with: bash scripts/setup_databases.sh" >&2
            exit 1
        fi
        skani search -d {params.sketch} \
                     -q {input.fnas} \
                     -o $tmp/raw.tsv \
                     -t {threads} --min-af 50

        python3 - <<PY
from pathlib import Path
from collections import defaultdict
import csv
rows = list(csv.DictReader(open("$tmp/raw.tsv"), delimiter="\t"))
nearest = defaultdict(float); n_close = defaultdict(int)
for r in rows:
    qf = Path(r.get("Query_file","")).stem
    if qf.endswith("_genomic"): qf = qf[:-len("_genomic")]
    try:
        ani = float(r.get("ANI", 0)); af  = float(r.get("Align_fraction_query", 0))
    except ValueError: continue
    if af < 50: continue
    if ani > nearest[qf]: nearest[qf] = ani
    if ani >= 95.0: n_close[qf] += 1
with open("{output.tsv}","w") as fh:
    fh.write("sample\tnearest_other_ani\tn_close_train\n")
    for q in sorted(nearest):
        fh.write(f"{{q}}\t{{nearest[q]:.3f}}\t{{n_close[q]}}\n")
PY
        """

# Keep GTDB-Tk available but optional
rule gtdbtk_classify:
    benchmark:
        OUT/"benchmarks"/"gtdbtk_classify.tsv"
    input:
        gd = OUT/"taxonomy"/"genome_dir",
    output:
        bac = OUT/"taxonomy"/"gtdbtk.bac120.summary.tsv",
    threads: config["threads"]["gtdbtk"]
    conda: "../envs/gtdbtk.yaml"
    params:
        db = config["gtdbtk_db"],
        out = OUT/"taxonomy",
    shell:
        r"""
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf --genome_dir {input.gd} --extension fna \
                           --out_dir {params.out} --cpus {threads} \
                           --pplacer_cpus 1 --skip_ani_screen
        touch {output.bac}
        """
