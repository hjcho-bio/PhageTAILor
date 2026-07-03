"""Build per-genome TSV + browser-friendly HTML report (v0.1.0).

Input:
  predictions.tsv  — per-region predictions from the classify rule.

Output:
  per_genome_predictions.tsv — one row per genome with max per-class scores.
  summary.html               — self-contained HTML viewable in any browser.
"""
import csv, html, sys
from datetime import datetime
from pathlib import Path
import pandas as pd

PRED = Path(snakemake.input.predictions)
OUT_TSV  = Path(snakemake.output.tsv)
OUT_HTML = Path(snakemake.output.html)
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

CLASSES = ["tailocin","Phage","T6SS","eCIS"]

df = pd.read_csv(PRED, sep="\t")
if df.empty:
    OUT_TSV.write_text("sample\tn_regions\t" + "\t".join(f"max_prob_{c}" for c in CLASSES) + "\n")
    OUT_HTML.write_text("<html><body><h1>PhageTAILor v0.1.0</h1><p>No candidate regions detected.</p></body></html>")
    sys.exit(0)

# Per-genome aggregates
agg_dict = {f"prob_{c}": "max" for c in CLASSES if f"prob_{c}" in df.columns}
for h in ["tailocin_binary","t6ss_binary","ecis_binary"]:
    if h in df.columns: agg_dict[h] = "max"
agg = df.groupby("sample").agg(agg_dict).reset_index()
agg["n_regions"] = df.groupby("sample").size().reindex(agg["sample"]).values

# Per-genome confidence band (use the band of the genome's highest-prob region)
if "confidence" in df.columns and "prob_tailocin" in df.columns:
    best = df.loc[df.groupby("sample")["prob_tailocin"].idxmax()]
    agg["confidence"] = agg["sample"].map(best.set_index("sample")["confidence"].to_dict()).fillna("no_match")

# Genome-level call source per class, chosen from clean-negative panel behaviour:
#  - tailocin, T6SS -> binary heads @0.5. The 5-way multiclass softmax under-calls
#    at 0.5 (mass split across classes -> ~0.58 tailocin recall); the binary heads
#    are well-separated (tailocin precision 0.81 on clean negatives, 0 FP on hard
#    same-clade and structural negatives; T6SS in-distribution AUPRC 0.91).
#  - eCIS -> MULTICLASS prob_eCIS. The eCIS binary head saturates at ~1.0 under its
#    extreme imbalance (scale_pos_weight >1000) and over-calls badly (8/9 structural
#    negatives falsely positive), whereas the multiclass head correctly assigns those
#    ~1e-8 eCIS probability.
#  - Phage -> multiclass (no dedicated binary head; prophage also covered by geNomad).
SCORE_SRC = {"tailocin":"tailocin_binary","T6SS":"t6ss_binary","eCIS":"prob_eCIS","Phage":"prob_Phage"}
for c in CLASSES:
    src = SCORE_SRC.get(c, f"prob_{c}")
    agg[f"score_{c}"] = agg[src] if src in agg.columns else agg.get(f"prob_{c}", 0.0)
    agg[f"call_{c}"] = (agg[f"score_{c}"] >= 0.5).astype(int)

agg = agg.sort_values("sample")
agg.to_csv(OUT_TSV, sep="\t", index=False)

# HTML
n = len(agg)
counts = {c: int(agg.get(f"call_{c}", pd.Series([0]*n)).sum()) for c in CLASSES}
bands = agg.get("confidence", pd.Series(["unknown"]*n)).value_counts().to_dict()

def badge(b):
    colors = {"high":"#2ecc71","medium":"#f1c40f","low":"#e67e22","no_match":"#95a5a6"}
    color = colors.get(b, "#bdc3c7")
    return f'<span style="background:{color};color:white;padding:2px 8px;border-radius:10px;font-size:0.8em">{html.escape(str(b))}</span>'

def cell(p):
    try: v = float(p)
    except: v = 0.0
    bg = "#2ecc71" if v >= 0.5 else ("#f1c40f" if v >= 0.3 else "#ecf0f1")
    fg = "white" if v >= 0.5 else "#2c3e50"
    return f'<td style="background:{bg};color:{fg};text-align:right;padding:4px 8px">{v:.3f}</td>'

parts = [f'''<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>PhageTAILor v0.1.0 report</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; max-width: 1200px; margin: 40px auto; padding: 0 20px; color: #2c3e50; }}
h1 {{ border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
h2 {{ color: #3498db; margin-top: 30px; }}
.grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; margin: 20px 0; }}
.box {{ background: #ecf0f1; padding: 15px; border-radius: 5px; }}
.box h3 {{ margin: 0 0 5px; color: #34495e; font-size: 0.9em; }}
.box .big {{ font-size: 1.8em; font-weight: bold; color: #2c3e50; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.9em; }}
th {{ background: #34495e; color: white; padding: 8px; text-align: left; }}
td {{ border-bottom: 1px solid #ecf0f1; padding: 6px 8px; }}
tr:hover {{ background: #f8f9fa; }}
.footer {{ margin-top: 40px; color: #7f8c8d; font-size: 0.85em; text-align: center; }}
.note {{ background: #fff7e6; border-left: 4px solid #f39c12; padding: 10px 15px; margin: 15px 0; font-size: 0.9em; }}
</style></head><body>

<h1>PhageTAILor v0.1.0 — results summary</h1>
<p>Generated {datetime.now().strftime("%Y-%m-%d %H:%M")} &middot; <b>{n}</b> genomes analyzed</p>

<h2>Per-PTE-class call counts (score &ge; 0.5; binary heads for tailocin/T6SS/eCIS, multiclass for Phage)</h2>
<div class="grid">
  <div class="box"><h3>tailocin</h3><div class="big">{counts.get("tailocin",0)}</div><div>{100*counts.get("tailocin",0)/max(n,1):.1f}% of genomes</div></div>
  <div class="box"><h3>Phage (prophage)</h3><div class="big">{counts.get("Phage",0)}</div><div>{100*counts.get("Phage",0)/max(n,1):.1f}%</div></div>
  <div class="box"><h3>T6SS</h3><div class="big">{counts.get("T6SS",0)}</div><div>{100*counts.get("T6SS",0)/max(n,1):.1f}%</div></div>
  <div class="box"><h3>eCIS</h3><div class="big">{counts.get("eCIS",0)}</div><div>{100*counts.get("eCIS",0)/max(n,1):.1f}%</div></div>
</div>

<h2>Confidence band distribution</h2>
<p>How close each query genome is to the v0.1.0 training pool. Predictions on <code>high</code>-confidence genomes are most trustworthy; <code>no_match</code> means the query is in a clade the model has never seen.</p>
<div class="grid">
  <div class="box"><h3>high (ANI &ge; 95)</h3><div class="big">{bands.get("high",0)}</div></div>
  <div class="box"><h3>medium (90&ndash;95)</h3><div class="big">{bands.get("medium",0)}</div></div>
  <div class="box"><h3>low (&lt; 90)</h3><div class="big">{bands.get("low",0)}</div></div>
  <div class="box"><h3>no_match</h3><div class="big">{bands.get("no_match",0)}</div></div>
</div>

<h2>Per-genome results (top {min(100,n)} by tailocin score)</h2>
<div class="note">Green = score &ge; 0.5 (positive call) &middot; Yellow = 0.3&ndash;0.5 (borderline) &middot; Gray = &lt; 0.3 (no signal)</div>

<table>
<thead><tr>
  <th>Sample</th><th>Confidence</th>
  <th style="text-align:right">tailocin</th>
  <th style="text-align:right">Phage</th>
  <th style="text-align:right">T6SS</th>
  <th style="text-align:right">eCIS</th>
  <th>n_regions</th>
</tr></thead><tbody>
''']

sort_key = "score_tailocin" if "score_tailocin" in agg.columns else "sample"
for _, r in agg.sort_values(sort_key, ascending=False).head(100).iterrows():
    parts.append(
        f"<tr><td><code>{html.escape(str(r['sample']))}</code></td><td>{badge(r.get('confidence','no_match'))}</td>"
        f"{cell(r.get('score_tailocin',0))}{cell(r.get('score_Phage',0))}"
        f"{cell(r.get('score_T6SS',0))}{cell(r.get('score_eCIS',0))}"
        f"<td style='text-align:right'>{int(r.get('n_regions',0))}</td></tr>"
    )

parts.append(f'''
</tbody></table>

<p style="margin-top:15px">
<b>For the full table</b>, open <code>per_genome_predictions.tsv</code> in Excel.<br>
<b>For per-region details</b>, open <code>../classify/predictions.tsv</code>.<br>
<b>For genome-browser visualization</b>, load <code>../features/regions.bed</code> in IGV / JBrowse / UCSC.
</p>

<div class="footer">
PhageTAILor v0.1.0 &middot; multiclass head: LightGBM 3.3 &middot; 5-detector union (geNomad + SecReT6 + eCIStem + tail-gene + PHROGS-tail)<br>
For step-by-step documentation, run <code>phagetailor explain &lt;step&gt;</code>.
</div>

</body></html>''')

OUT_HTML.write_text("".join(parts))
print(f"[report] wrote {OUT_TSV} and {OUT_HTML} ({n} genomes)", file=sys.stderr)
