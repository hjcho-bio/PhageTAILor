"""Apply the LightGBM models (multiclass + 3 binary heads) to per-region features,
with isotonic calibration of the binary heads.

Loads from `resources/model/`:
  multiclass.txt           (5-way: Phage, tailocin, T6SS, eCIS, negative)
  tailocin_binary.txt
  t6ss_binary.txt
  ecis_binary.txt
  feature_columns.json     (column order shared by all heads)
  binary_calibrators.pkl   (per-head IsotonicRegression; applied per-region
                            before per-genome aggregation)

The binary-head outputs reported downstream (`tailocin_binary`, `t6ss_binary`,
`ecis_binary`) are the CALIBRATED probabilities; the raw uncalibrated outputs
are retained as `*_raw` for diagnostics.
"""
import json, pickle, sys
from pathlib import Path
import numpy as np
import pandas as pd
import lightgbm as lgb

FEATS = Path(snakemake.input.features) if hasattr(snakemake.input,"features") else Path(snakemake.input[0])
MODEL_DIR = Path(snakemake.config["model_dir"])
OUT = Path(snakemake.output.predictions)
OUT.parent.mkdir(parents=True, exist_ok=True)

with open(MODEL_DIR/"feature_columns.json") as fh:
    FEAT_COLS = json.load(fh)

CLASSES = ["negative","tailocin","Phage","T6SS","eCIS"]

def load_booster(name):
    return lgb.Booster(model_file=str(MODEL_DIR/f"{name}.txt"))

m_mc       = load_booster("multiclass")
m_tailocin = load_booster("tailocin_binary")
m_t6ss     = load_booster("t6ss_binary")
m_ecis     = load_booster("ecis_binary")

# Per-head isotonic calibrators (fit on held-out 5-fold CV predictions).
CAL_PATH = MODEL_DIR/"binary_calibrators.pkl"
if CAL_PATH.exists():
    with open(CAL_PATH, "rb") as fh:
        CALIBRATORS = pickle.load(fh)
    print(f"[classify] loaded isotonic calibrators for {list(CALIBRATORS.keys())}", file=sys.stderr)
else:
    CALIBRATORS = {}
    print("[classify] WARNING: binary_calibrators.pkl not found; reporting raw binary-head scores", file=sys.stderr)

# Load features
try:
    df = pd.read_parquet(FEATS)
except Exception:
    alt = str(FEATS).replace(".parquet",".tsv.gz")
    df = pd.read_csv(alt, sep="\t", compression="infer", low_memory=False)

# Align to model schema (zero-fill missing features)
for c in FEAT_COLS:
    if c not in df.columns: df[c] = 0.0
X = df[FEAT_COLS].astype(np.float32).values

# Multiclass softmax (LightGBM returns shape [n, n_class])
proba_mc = m_mc.predict(X)
if proba_mc.ndim == 1:
    proba_mc = proba_mc.reshape(-1, len(CLASSES))
for i, c in enumerate(CLASSES):
    df[f"prob_{c}"] = proba_mc[:, i]
df["pred_class"] = [CLASSES[i] for i in proba_mc.argmax(axis=1)]
df["pred_score"] = proba_mc.max(axis=1)

# Binary heads — raw + isotonically-calibrated
def predict_binary(name, booster):
    raw = booster.predict(X)
    df[f"{name}_raw"] = raw
    iso = CALIBRATORS.get(name)
    df[name] = iso.predict(raw) if iso is not None else raw

predict_binary("tailocin_binary", m_tailocin)
predict_binary("t6ss_binary",     m_t6ss)
predict_binary("ecis_binary",     m_ecis)

# Confidence band
def confidence_band(ani):
    if ani is None or ani <= 0: return "no_match"
    if ani >= 95: return "high"
    if ani >= 90: return "medium"
    return "low"
df["confidence"] = df.get("nearest_other_ani", 0).apply(confidence_band)

# Output columns (calibrated binary heads as the reported scores; raw kept for diagnostics)
keep = ["sample","scaffold","start","end","src","candidate_id",
        "pred_class","pred_score",
        "prob_negative","prob_tailocin","prob_Phage","prob_T6SS","prob_eCIS",
        "tailocin_binary","t6ss_binary","ecis_binary",
        "tailocin_binary_raw","t6ss_binary_raw","ecis_binary_raw",
        "nearest_other_ani","n_close_train","confidence"]
df_out = df[[c for c in keep if c in df.columns]]
df_out.to_csv(OUT, sep="\t", index=False)
print(f"[classify] wrote {len(df_out)} predictions to {OUT}", file=sys.stderr)
