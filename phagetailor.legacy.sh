#!/bin/bash
# Thin CLI wrapper: phagetailor run --input <list.txt> --out <dir> [--modules ...]
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 run [options]

Required:
  --input <file>      Path to a text file with one .fna path per line
  --out   <dir>       Output directory

Optional:
  --modules <csv>     Comma-separated subset of modules (default: all)
                      Available: annotation,taxonomy,phage,tailocin,t6ss,ecis,
                                 features,classify,report
  --cores <int>       Number of CPU cores (default: 8)
  --use-conda         Use snakemake's conda env management
  --                  Pass remaining args directly to snakemake

Examples:
  $0 run --input genomes.txt --out results/ --cores 32
  $0 run --input genomes.txt --out results/ --modules tailocin,t6ss --cores 16
EOF
}

if [ $# -eq 0 ] || [ "$1" != "run" ]; then usage; exit 1; fi
shift

INPUT=""; OUT=""; MODULES="all"; CORES=8; USE_CONDA=""
EXTRA_ARGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --input)   INPUT="$2"; shift 2 ;;
    --out)     OUT="$2"; shift 2 ;;
    --modules) MODULES="$2"; shift 2 ;;
    --cores)   CORES="$2"; shift 2 ;;
    --use-conda) USE_CONDA="--use-conda"; shift ;;
    --) shift; EXTRA_ARGS=("$@"); break ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1"; usage; exit 1 ;;
  esac
done

[ -z "$INPUT" ] && { echo "ERROR: --input required"; exit 1; }
[ -z "$OUT" ]   && { echo "ERROR: --out required"; exit 1; }

PHAGETAILOR_DIR="$(cd "$(dirname "$0")" && pwd)"
exec snakemake -s "${PHAGETAILOR_DIR}/workflow/Snakefile" \
               --cores "${CORES}" \
               ${USE_CONDA} \
               --config "genome_list=${INPUT}" "output_dir=${OUT}" "modules=${MODULES}" \
               "${EXTRA_ARGS[@]}"
