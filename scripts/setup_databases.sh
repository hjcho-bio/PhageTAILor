#!/usr/bin/env bash
# Download the PhageTAILor reference DBs that are too large to store in Git.
# All are hosted on a single Zenodo record:
#   - tail-HMM library   -> resources/hmm_tail_db/                  (110 MB download)
#   - PHROGS mmseqs2 DB  -> resources/phrogs_db/phrogs_mmseqs_db/   ( 85 MB download)
#   - training sketch    -> resources/training_sketch/sketch_db/    (2.9 GB, split into
#                                                                    6 parts; confidence band)
#
# (The geNomad DB is third-party -- see setup_genomad_db.sh.)
#
# Usage:
#   bash scripts/setup_databases.sh            # get all three
#   bash scripts/setup_databases.sh sketch     # just the training sketch
#   bash scripts/setup_databases.sh reference  # just tail-HMM + PHROGS
set -euo pipefail

# Zenodo record 21152308 (https://doi.org/10.5281/zenodo.21152308).
# Override with $PHAGETAILOR_ZENODO_BASE if you mirror the tarballs elsewhere.
ZENODO_BASE="${PHAGETAILOR_ZENODO_BASE:-https://zenodo.org/records/21152308/files}"

WHAT="${1:-all}"
DIR=$(cd "$(dirname "$0")"/.. && pwd)

CURL_OPTS=(-fL --retry 5 --retry-delay 5 --retry-all-errors)

# md5 checksums of the published artefacts. Keep in sync with the Zenodo record.
md5_for() {
    case "$1" in
      phagetailor_hmm_tail_db_v0.1.0.tar.gz)    echo c196e2c966dec6271e4357cc76fa31a8 ;;
      phagetailor_phrogs_db_v0.1.0.tar.gz)      echo 92ac74caa6ceaea8aa1661e0561a9cdb ;;
      phagetailor_sketch_v0.1.0.tar.gz.part-aa) echo dd81466d71a2bc893be907c14f9266b2 ;;
      phagetailor_sketch_v0.1.0.tar.gz.part-ab) echo 3fe0d995ef3d3bdd5a08a2fc486bd089 ;;
      phagetailor_sketch_v0.1.0.tar.gz.part-ac) echo 4f8e98c4bcadb4393081ec50bd14c4b9 ;;
      phagetailor_sketch_v0.1.0.tar.gz.part-ad) echo 56f9b779879171fee99e44ea2f2d35d3 ;;
      phagetailor_sketch_v0.1.0.tar.gz.part-ae) echo dfa112b49179851e5d14d7167f8007dd ;;
      phagetailor_sketch_v0.1.0.tar.gz.part-af) echo e7bb71dec6797d038f5cb931089bf3d8 ;;
      *) echo "" ;;
    esac
}

md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
    else md5 -q "$1"; fi
}

verify() {
    local path="$1" name="$2" want; want=$(md5_for "$name")
    [ -z "$want" ] && return 0
    local got; got=$(md5_of "$path")
    if [ "$got" != "$want" ]; then
        echo "error: checksum mismatch for ${name}" >&2
        echo "       expected ${want}" >&2
        echo "       got      ${got}" >&2
        echo "       The download was truncated or corrupted. Delete it and re-run." >&2
        exit 1
    fi
}

download() {   # download() URL DEST -- with resume
    local name="$1" dest="$2"
    echo "Downloading ${name} ..."
    curl "${CURL_OPTS[@]}" -C - "${ZENODO_BASE}/${name}" -o "$dest"
    verify "$dest" "$name"
}

fetch() {      # single-file tarball
    local file="$1" dest="$2" sentinel="$3"
    if [ -e "${DIR}/${sentinel}" ]; then
        echo "[skip] ${sentinel} already present."
        return
    fi
    local tmp; tmp=$(mktemp -d); trap 'rm -rf "$tmp"' RETURN
    download "$file" "$tmp/a.tar.gz"
    echo "Extracting into ${DIR}/${dest} ..."
    mkdir -p "${DIR}/${dest}"
    tar -xzf "$tmp/a.tar.gz" -C "${DIR}/${dest}"
    [ -e "${DIR}/${sentinel}" ] || { echo "error: expected ${sentinel} after extraction" >&2; exit 1; }
}

fetch_split() {   # multi-part tarball, reassembled on the fly
    local base="$1" dest="$2" sentinel="$3"; shift 3
    local parts=("$@")
    if [ -e "${DIR}/${sentinel}" ]; then
        echo "[skip] ${sentinel} already present."
        return
    fi
    local tmp; tmp=$(mktemp -d); trap 'rm -rf "$tmp"' RETURN
    echo "The training sketch is ~2.9 GB, downloaded as ${#parts[@]} parts."
    local n=0
    for suf in "${parts[@]}"; do
        n=$((n+1))
        echo "  [part ${n}/${#parts[@]}]"
        download "${base}.part-${suf}" "$tmp/${suf}"
    done
    echo "Reassembling and extracting into ${DIR}/${dest} ..."
    mkdir -p "${DIR}/${dest}"
    local ordered=()
    for suf in "${parts[@]}"; do ordered+=("$tmp/${suf}"); done
    cat "${ordered[@]}" | tar -xz -C "${DIR}/${dest}"
    [ -e "${DIR}/${sentinel}" ] || { echo "error: expected ${sentinel} after extraction" >&2; exit 1; }
}

if [ "$WHAT" = "all" ] || [ "$WHAT" = "reference" ]; then
    fetch phagetailor_hmm_tail_db_v0.1.0.tar.gz "resources"           "resources/hmm_tail_db/hmm_tail_full.hmm"
    fetch phagetailor_phrogs_db_v0.1.0.tar.gz   "resources/phrogs_db" "resources/phrogs_db/phrogs_mmseqs_db/phrogs_profile_db_seq"
fi
if [ "$WHAT" = "all" ] || [ "$WHAT" = "sketch" ]; then
    fetch_split phagetailor_sketch_v0.1.0.tar.gz "resources/training_sketch" \
                "resources/training_sketch/sketch_db" aa ab ac ad ae af
fi

echo "Done."
