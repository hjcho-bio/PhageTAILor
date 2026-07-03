#!/usr/bin/env bash
# Download the PhageTAILor reference DBs that are too large to store in Git.
# All are hosted on a single Zenodo record:
#   - tail-HMM library   -> resources/hmm_tail_db/                  (~350 MB)
#   - PHROGS mmseqs2 DB   -> resources/phrogs_db/phrogs_mmseqs_db/  (~400 MB)
#   - training sketch     -> resources/training_sketch/sketch_db/   (~2.7 GB, confidence band)
#
# (The geNomad DB is third-party — see setup_genomad_db.sh.)
#
# Usage:
#   bash scripts/setup_databases.sh            # Get all three
#   bash scripts/setup_databases.sh sketch     # Just the training sketch
#   bash scripts/setup_databases.sh reference  # Just tail-HMM + PHROGS
set -euo pipefail

# Zenodo record 21152308 (https://doi.org/10.5281/zenodo.21152308).
# Override with $PHAGETAILOR_ZENODO_BASE if you mirror the tarballs elsewhere.
ZENODO_BASE="${PHAGETAILOR_ZENODO_BASE:-https://zenodo.org/records/21152308/files}"

WHAT="${1:-all}"
DIR=$(cd "$(dirname "$0")"/.. && pwd)

if [ -z "$ZENODO_BASE" ]; then
    echo "error: no Zenodo base URL. Set it in this script or via PHAGETAILOR_ZENODO_BASE:" >&2
    echo "  PHAGETAILOR_ZENODO_BASE=https://zenodo.org/records/<ID>/files bash scripts/setup_databases.sh" >&2
    exit 1
fi

fetch() {
    local file="$1" dest="$2" sentinel="$3"
    if [ -e "${DIR}/${sentinel}" ]; then
        echo "[skip] ${sentinel} already present."
        return
    fi
    local tmp; tmp=$(mktemp -d); trap 'rm -rf "$tmp"' RETURN
    echo "Downloading ${file} ..."
    if command -v curl >/dev/null 2>&1; then curl -fL "${ZENODO_BASE}/${file}" -o "$tmp/a.tar.gz"
    else wget -O "$tmp/a.tar.gz" "${ZENODO_BASE}/${file}"; fi
    echo "Extracting into ${DIR}/${dest} ..."
    mkdir -p "${DIR}/${dest}"
    tar -xzf "$tmp/a.tar.gz" -C "${DIR}/${dest}"
    [ -e "${DIR}/${sentinel}" ] || { echo "error: expected ${sentinel} after extraction" >&2; exit 1; }
}

if [ "$WHAT" = "all" ] || [ "$WHAT" = "reference" ]; then
    fetch phagetailor_hmm_tail_db_v0.1.0.tar.gz "resources"           "resources/hmm_tail_db/hmm_tail_full.hmm"
    fetch phagetailor_phrogs_db_v0.1.0.tar.gz   "resources/phrogs_db" "resources/phrogs_db/phrogs_mmseqs_db/phrogs_profile_db_seq"
fi
if [ "$WHAT" = "all" ] || [ "$WHAT" = "sketch" ]; then
    fetch phagetailor_sketch_v0.1.0.tar.gz "resources/training_sketch" "resources/training_sketch/sketch_db"
fi

echo "Done."