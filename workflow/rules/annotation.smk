"""Module: annotation — Pyrodigal-gv proteins + MMseqs2 vs four sequence DBs + HMMER vs the tail-HMM library.

Pyrodigal-gv is the same gene caller used at training time (single-tool consistency
across train and inference datasets). The four MMseqs2 searches feed the cluster-based
detectors in `tailocin`, `t6ss`, and `ecis`, plus the PHROGS category counts in `features`.
The HMMER search feeds the HMM-tail cluster detector (`rules/hmm_tail.smk`) and supplies
the region- and genome-level HMM features in `features`.
"""

rule pyrodigal_call_proteins:
    benchmark:
        OUT/"benchmarks"/"pyrodigal_call_proteins.{s}.tsv"
    input:
        fna = lambda wc: GENOMES[wc.s],
    output:
        faa = OUT/"annotation"/"{s}.faa",
        gff = OUT/"annotation"/"{s}.gff",
    threads: config["threads"]["pyrodigal"]
    conda: "../envs/annotation.yaml"
    shell:
        r"""
        pyrodigal-gv -p meta -j {threads} -i {input.fna} -a {output.faa} -f gff -o {output.gff}
        """

# Generic MMseqs2 search: written once, instantiated four times below.
def _mmseqs_search_shell(ref_key, max_seqs, sensitivity):
    return r"""
        tmp=$(mktemp -d)
        mmseqs createdb {input.query} $tmp/q --shuffle 0
        mmseqs search $tmp/q {params.ref} $tmp/r $tmp/tmp \
               -e {params.evalue} --threads {threads} --max-seqs %d -s %s
        mmseqs convertalis $tmp/q {params.ref} $tmp/r {output.m8} \
               --format-output "query,target,pident,alnlen,evalue,bits" \
               --threads {threads}
        rm -rf $tmp
        """ % (max_seqs, sensitivity)

rule mmseqs_search_secret6:
    benchmark:
        OUT/"benchmarks"/"mmseqs_search_secret6.{s}.tsv"
    input:  query = OUT/"annotation"/"{s}.faa"
    output: m8 = OUT/"annotation"/"{s}.secret6.m8"
    threads: config["threads"]["mmseqs"]
    conda: "../envs/annotation.yaml"
    params: ref = config["secret6_mmseqs_db"], evalue = config["mmseqs_evalue"]
    shell:  _mmseqs_search_shell("secret6_mmseqs_db", 5, "5.7")

rule mmseqs_search_ecistem:
    benchmark:
        OUT/"benchmarks"/"mmseqs_search_ecistem.{s}.tsv"
    input:  query = OUT/"annotation"/"{s}.faa"
    output: m8 = OUT/"annotation"/"{s}.ecistem.m8"
    threads: config["threads"]["mmseqs"]
    conda: "../envs/annotation.yaml"
    params: ref = config["ecistem_mmseqs_db"], evalue = config["mmseqs_evalue"]
    shell:  _mmseqs_search_shell("ecistem_mmseqs_db", 5, "5.7")

rule mmseqs_search_phrogs:
    benchmark:
        OUT/"benchmarks"/"mmseqs_search_phrogs.{s}.tsv"
    input:  query = OUT/"annotation"/"{s}.faa"
    output: m8 = OUT/"annotation"/"{s}.phrogs.m8"
    threads: config["threads"]["mmseqs"]
    conda: "../envs/annotation.yaml"
    params: ref = config["phrogs_mmseqs_db"], evalue = config["mmseqs_evalue"]
    shell:  _mmseqs_search_shell("phrogs_mmseqs_db", 3, "5")

rule mmseqs_search_tailgenes:
    benchmark:
        OUT/"benchmarks"/"mmseqs_search_tailgenes.{s}.tsv"
    input:  query = OUT/"annotation"/"{s}.faa"
    output: m8 = OUT/"annotation"/"{s}.tailgenes.m8"
    threads: config["threads"]["mmseqs"]
    conda: "../envs/annotation.yaml"
    params: ref = config["tailgenes_mmseqs_db"], evalue = config["mmseqs_evalue"]
    shell:  _mmseqs_search_shell("tailgenes_mmseqs_db", 5, "5.7")

rule hmmsearch_tail:
    benchmark:
        OUT/"benchmarks"/"hmmsearch_tail.{s}.tsv"
    """HMMER search against the bundled tail-HMM library. The output .tbl file
    feeds the HMM-tail clusterer (rules/hmm_tail.smk) and the HMM features (features.smk).
    The full 597-profile library is used so the per-region feature builder can read
    capsid/portal hits as discriminative module-count features, while the cluster
    detector internally restricts to the 312-profile cassette-only subset.
    """
    input:  query = OUT/"annotation"/"{s}.faa"
    output: tbl = OUT/"annotation"/"{s}.hmmtail.tbl"
    threads: config["threads"]["hmmer"]
    conda: "../envs/annotation.yaml"
    params:
        db     = config["hmm_tail_db"],
        evalue = config["hmmer_evalue"],
    shell:
        r"""
        hmmsearch --cpu {threads} -E {params.evalue} --noali \
                  --tblout {output.tbl} {params.db} {input.query} >/dev/null
        """
