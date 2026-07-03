#!/usr/bin/env bash
# Download the geNomad reference database (~5 GB, one-time).
#
# Usage:
#   bash scripts/setup_genomad_db.sh [DEST_DIR]
#
# DEST_DIR defaults to ./resources. The database is written to
# DEST_DIR/genomad_db. After it finishes, point `genomad_db` in
# config/config.yaml at that path.
set -euo pipefail

DEST_DIR="${1:-resources}"

if ! command -v genomad >/dev/null 2>&1; then
    echo "error: 'genomad' not found on PATH." >&2
    echo "Activate the phagetailor conda env first: conda activate phagetailor" >&2
    exit 1
fi

mkdir -p "$DEST_DIR"
echo "Downloading geNomad database into ${DEST_DIR}/genomad_db ..."
genomad download-database "$DEST_DIR"

echo
echo "Done. Set this in config/config.yaml:"
echo "  genomad_db: \"$(cd "$DEST_DIR" && pwd)/genomad_db\""