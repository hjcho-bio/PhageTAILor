#!/usr/bin/env bash
# Download the geNomad reference database (~5 GB, one-time).
#
# This does NOT require geNomad to be installed. The database is a plain
# tarball published by the geNomad authors on Zenodo; we fetch it directly.
# If geNomad happens to be on PATH we use it, since it handles mirroring.
#
# Usage:
#   bash scripts/setup_genomad_db.sh [DEST_DIR]
#
# DEST_DIR defaults to ./resources. The database is written to
# DEST_DIR/genomad_db, and config/config.yaml is updated to point at it.
set -euo pipefail

DEST_DIR="${1:-resources}"
DIR=$(cd "$(dirname "$0")"/.. && pwd)
GENOMAD_DB_URL="${PHAGETAILOR_GENOMAD_DB_URL:-https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz}"

mkdir -p "$DEST_DIR"
ABS_DEST=$(cd "$DEST_DIR" && pwd)

if [ -e "${ABS_DEST}/genomad_db/genomad_db" ] || [ -d "${ABS_DEST}/genomad_db" ]; then
    echo "[skip] ${ABS_DEST}/genomad_db already present."
else
    if command -v genomad >/dev/null 2>&1; then
        echo "geNomad found on PATH; using 'genomad download-database'."
        genomad download-database "$ABS_DEST"
    else
        echo "geNomad not on PATH -- downloading the database tarball directly."
        echo "  ${GENOMAD_DB_URL}"
        tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
        curl -fL --retry 5 --retry-delay 5 --retry-all-errors -C - \
             "$GENOMAD_DB_URL" -o "$tmp/genomad_db.tar.gz"
        echo "Extracting into ${ABS_DEST} ..."
        tar -xzf "$tmp/genomad_db.tar.gz" -C "$ABS_DEST"
    fi
fi

if [ ! -d "${ABS_DEST}/genomad_db" ]; then
    echo "error: expected ${ABS_DEST}/genomad_db after download." >&2
    exit 1
fi

# Point config/config.yaml at the database we just installed.
CFG="${DIR}/config/config.yaml"
if [ -f "$CFG" ]; then
    python3 - "$CFG" "${ABS_DEST}/genomad_db" <<'PY'
import re, sys
cfg, path = sys.argv[1], sys.argv[2]
text = open(cfg).read()
new, n = re.subn(r'^genomad_db:.*$', f'genomad_db:  "{path}"', text, flags=re.M)
if n:
    open(cfg, 'w').write(new)
    print(f"Updated genomad_db in {cfg}")
else:
    print(f"Could not find a genomad_db line in {cfg}; set it manually to {path}")
PY
fi

echo
echo "Done. geNomad DB installed at ${ABS_DEST}/genomad_db"
