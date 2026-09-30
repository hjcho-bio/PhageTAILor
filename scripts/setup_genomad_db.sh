#!/usr/bin/env bash
# Download the geNomad reference database (~1.5 GB download, ~5 GB extracted).
#
# This does NOT require geNomad to be installed first. If geNomad is not on
# PATH we build a throwaway conda env for it and use its own downloader, so
# the database version always matches the geNomad release the pipeline pins.
# (Hardcoding a Zenodo URL here would silently rot when geNomad bumps its DB.)
#
# Usage:
#   bash scripts/setup_genomad_db.sh [DEST_DIR]
#
# DEST_DIR defaults to ./resources. The database is written to
# DEST_DIR/genomad_db, and config/config.yaml is pointed at it.
set -euo pipefail

DEST_DIR="${1:-resources}"
DIR=$(cd "$(dirname "$0")"/.. && pwd)
ENV_FILE="${DIR}/workflow/envs/phage.yaml"

mkdir -p "$DEST_DIR"
ABS_DEST=$(cd "$DEST_DIR" && pwd)

if [ -d "${ABS_DEST}/genomad_db" ]; then
    echo "[skip] ${ABS_DEST}/genomad_db already present."
elif command -v genomad >/dev/null 2>&1; then
    echo "geNomad found on PATH; using 'genomad download-database'."
    genomad download-database "$ABS_DEST"
else
    CONDA_EXE_BIN="$(command -v mamba || command -v conda || true)"
    if [ -z "$CONDA_EXE_BIN" ]; then
        echo "error: neither geNomad nor conda/mamba is available." >&2
        echo "  Install conda (miniforge), or install geNomad yourself and re-run." >&2
        exit 1
    fi
    TMP_ENV="${DIR}/.genomad_dl_env"
    if [ ! -x "${TMP_ENV}/bin/genomad" ]; then
        echo "Building a temporary env to obtain geNomad (from ${ENV_FILE}) ..."
        "$CONDA_EXE_BIN" env create -p "$TMP_ENV" -f "$ENV_FILE" -y >/dev/null
    fi
    echo "Downloading the geNomad database ..."
    "${TMP_ENV}/bin/genomad" download-database "$ABS_DEST"
    echo "Removing the temporary env ..."
    "$CONDA_EXE_BIN" env remove -p "$TMP_ENV" -y >/dev/null 2>&1 || rm -rf "$TMP_ENV"
fi

if [ ! -d "${ABS_DEST}/genomad_db" ]; then
    echo "error: expected ${ABS_DEST}/genomad_db after download." >&2
    exit 1
fi

CFG="${DIR}/config/config.yaml"
if [ -f "$CFG" ]; then
    python3 - "$CFG" "${ABS_DEST}/genomad_db" <<'PY'
import re, sys
cfg, path = sys.argv[1], sys.argv[2]
text = open(cfg).read()
new, n = re.subn(r'^genomad_db:.*$', f'genomad_db:  "{path}"', text, flags=re.M)
if n:
    open(cfg, 'w').write(new); print(f"Updated genomad_db in {cfg}")
else:
    print(f"Could not find a genomad_db line in {cfg}; set it manually to {path}")
PY
fi

echo
echo "Done. geNomad DB installed at ${ABS_DEST}/genomad_db"
