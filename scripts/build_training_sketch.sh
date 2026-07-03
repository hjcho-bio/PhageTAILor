#!/bin/bash
# Rebuild the skani training sketch from a local directory of training genomes.
#
# Most users do NOT need this — download the prebuilt sketch instead:
#   bash scripts/setup_databases.sh
# This script is for reproducibility: it re-sketches the training genomes named
# in resources/training_sketch/training_accessions.txt (5,989 NCBI GCF/GCA +
# 512 IMG/JGI genomes). You must obtain those FASTA files yourself — NCBI ones
# via `datasets download genome accession`, IMG ones from https://img.jgi.doe.gov/ —
# and place them in GENOME_DIR (filenames must start with the accession).
#
# Usage: scripts/build_training_sketch.sh GENOME_DIR [THREADS]
set -euo pipefail

GENOME_DIR=${1:?usage: build_training_sketch.sh GENOME_DIR [THREADS]}
THREADS=${2:-32}
DIR=$(cd "$(dirname "$0")"/.. && pwd)
MANIFEST=${DIR}/resources/training_sketch/training_accessions.txt
OUT=${DIR}/resources/training_sketch/sketch_db

[ -f "${MANIFEST}" ] || { echo "missing ${MANIFEST}"; exit 1; }
[ -d "${GENOME_DIR}" ] || { echo "genome dir not found: ${GENOME_DIR}"; exit 1; }
[ -e "${OUT}" ] && { echo "sketch already exists at ${OUT}; delete it first to rebuild"; exit 0; }
command -v skani >/dev/null 2>&1 || { echo "skani not on PATH (e.g. conda activate gtdbtk)"; exit 1; }

# Resolve each accession to a FASTA file in GENOME_DIR (filename must start with the accession).
list=$(mktemp); trap 'rm -f "$list"' EXIT
missing=0
while IFS=$'\t' read -r acc _src; do
    [ "$acc" = "accession" ] && continue          # skip header
    f=$(find "${GENOME_DIR}" -maxdepth 1 \( -name "${acc}*.fna" -o -name "${acc}*.fasta" \) | head -1)
    if [ -n "$f" ]; then echo "$f" >> "$list"; else echo "  [missing] ${acc}" >&2; missing=$((missing+1)); fi
done < "${MANIFEST}"

n=$(wc -l < "$list")
echo "[$(date)] resolved ${n} genomes (${missing} missing) -> sketching to ${OUT}"
[ "$n" -gt 0 ] || { echo "no genomes resolved; check GENOME_DIR"; exit 1; }
skani sketch -l "$list" -o "${OUT}" -t "${THREADS}"
echo "[$(date)] done"