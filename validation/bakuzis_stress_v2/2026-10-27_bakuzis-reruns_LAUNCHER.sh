#!/usr/bin/env bash
## 2026-10-27_bakuzis-reruns_LAUNCHER.sh
## Turnkey launcher for the three reruns recommended in 20261027_fvs-modern_fvs-conus_bakuzis_assessment.md.
## Runs on OSC Cardinal (engine harness and projector metrics live there). Each phase is opt-in and gated;
## nothing is deleted and nothing is pushed. Phase outputs land under $OUT.
##
##   bash 2026-10-27_bakuzis-reruns_LAUNCHER.sh            # print the plan, check paths, run nothing
##   PHASE1=1 bash ...   # engine harness v2, GOMPIT and EXP-HAZ libs, both scored
##   PHASE2=1 bash ...   # FIA-stand ensemble seeding audit (year-0 BA vs label) on the banked ensemble + treeinit writer smoke
##   PHASE3=1 bash ...   # projector v4 metrics rescored with FIA percentiles and SDImax
set -uo pipefail
W=${W:-/fs/scratch/PUOM0008/crsfaaron/wt-engine}
ROOT=${ROOT:-/users/PUOM0008/crsfaaron/fvs-conus}
SCR=${SCR:-/fs/scratch/PUOM0008/crsfaaron/fvs_stress}
OUT=${OUT:-$SCR/bakuzis_v2_20261027}
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT"
source /etc/profile.d/lmod.sh 2>/dev/null; module load gcc/12.3.0 R/4.5.2 2>/dev/null || module load gcc/12.3.0 R/4.4.0

echo "== plan =="
echo "  scorer:   $HERE/bakuzis_stress_v2.R"
echo "  harness:  $HERE/engine_bakuzis_v2.py"
echo "  OUT:      $OUT"
for f in "$HERE/bakuzis_stress_v2.R" "$HERE/engine_bakuzis_v2.py" "$HERE/fia_ba_sdi_percentiles.R"; do [ -f "$f" ] && echo "  ok  $f" || { echo "  MISSING $f"; exit 1; }; done
LIB_G=${LIB_G:-$W/lib_gompit/FVSne.so}; LIB_E=${LIB_E:-$W/lib_exp/FVSne.so}
for f in "$LIB_G" "$LIB_E"; do [ -f "$f" ] && echo "  ok  $f" || echo "  (phase 1 needs) MISSING $f"; done
PAIRS=${PAIRS:-$ROOT/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds}
[ -f "$PAIRS" ] && echo "  ok  $PAIRS" || echo "  (phase 3 needs) MISSING $PAIRS"
ls "$SCR"/conus_eq_proj/*_metrics.csv 2>/dev/null | head -3 | sed 's/^/  metrics: /' || echo "  (phase 3 needs) no *_metrics.csv under $SCR/conus_eq_proj"

if [ "${PHASE1:-0}" = "1" ]; then
  echo "== phase 1: engine harness v2 (SITECODE, seeding gate), two survival libs =="
  export FVS_PROJECT_ROOT="$HOME/fvs-modern" FVS_CONFIG_DIR="$W/config" FVS_GOMPIT=1 FVS_GOMPIT_COEF="$W/config/greg_mortality_coefficients.csv"
  export PYTHONPATH="$W:$FVS_PROJECT_ROOT:$FVS_PROJECT_ROOT/deployment/fvs2py:${PYTHONPATH:-}" PYTHONNOUSERSITE=1
  "$W/fvsvenv/bin/python" "$HERE/engine_bakuzis_v2.py" --wt "$W" --lib "$LIB_G" --tag gompit --lib "$LIB_E" --tag exp --out "$OUT/engine_v2" 2>&1 | tee "$OUT/engine_v2.log"
  grep -E 'GATE SUMMARY|BAKUZIS_ENGINE_V2_DONE' "$OUT/engine_v2.log" || { echo "phase 1 did not finish"; exit 1; }
  printf 'group,sdimax\nSpruce-Fir,1250\nNorthern-Hardwood,1384\nPine,1307\nOak-Pine,1384\n' > "$OUT/sdimax_ne_engine_tph.csv"   # engine SDIMAX table values x 2.471
  for t in gompit exp; do
    f="$OUT/engine_v2/engine_bakuzis_100yr_$t.csv"; [ -f "$f" ] || continue
    Rscript "$HERE/bakuzis_stress_v2.R" --input="$f" --source=engine --out="$OUT/score_engine_$t" --sdimax="$OUT/sdimax_ne_engine_tph.csv" --label="engine_v2_$t" > "$OUT/score_engine_$t.log" 2>&1
    grep -E '^Reineke axis|Wrote outputs|Error' "$OUT/score_engine_$t.log"
  done
  echo "phase 1 gate: scenarios failing seed/site gates per lib:"; for t in gompit exp; do [ -f "$OUT/engine_v2/seed_site_gate_$t.csv" ] && awk -F, -v t=$t 'NR>1 && $NF!="OK"{n++} END{print "  "t": "n+0" not OK"}' "$OUT/engine_v2/seed_site_gate_$t.csv"; done
fi

if [ "${PHASE2:-0}" = "1" ]; then
  echo "== phase 2: FIA-stand ensemble seeding audit =="
  ENS=${ENS:-$HOME/fvs-modern/calibration/output/bakuzis/bakuzis_uncertainty_summary_long.csv}
  Rscript "$HERE/bakuzis_stress_v2.R" --input="$ENS" --source=uncert --out="$OUT/score_fiaensemble" --sdimax="$OUT/sdimax_ne_engine_tph.csv" --label="fiaensemble" > "$OUT/score_fiaensemble.log" 2>&1
  Rscript -e 'suppressMessages(library(data.table)); d<-fread(commandArgs(TRUE)[1]); lab<-c(Low=60,Medium=120,High=180)*0.2295684
    s<-d[age==min(age),.(ba0=median(BAPH),tph0=median(TPH)),by=.(arm,density)][,label:=lab[density]][,ratio:=round(ba0/label,2)]
    print(s[order(arm,density)]); cat("\nA ratio far below 1 means the TREEINIT records reached the engine without their expansion factors.\n")' "$OUT/score_fiaensemble/fiaensemble_long_metric.csv"
  echo "next: inspect calibration/python/inventory_keyfile.py and fia_stand_generator.py for the per-record TPA (expansion) written to TREEINIT; then re-run bakuzis_uncertainty_comparison.py with a year-0 gate."
fi

if [ "${PHASE3:-0}" = "1" ]; then
  echo "== phase 3: projector v4 metrics rescored with observed percentiles =="
  Rscript "$HERE/fia_ba_sdi_percentiles.R" "$PAIRS" "$OUT/fia_ba_sdi_percentiles.csv" || { echo "percentile script stopped; fix column names and rerun"; exit 1; }
  for f in "$SCR"/conus_eq_proj/*_metrics.csv; do
    v=$(basename "$f" _metrics.csv)
    var=$(echo "$v" | sed -E 's/^conus_eq_([a-z]+)_.*/\1/' | tr a-z A-Z)
    p999=$(awk -F, -v v="$var" 'NR>1 && toupper($1)==v {print $4}' "$OUT/fia_ba_sdi_percentiles.csv")
    sdi=$(awk -F, -v v="$var" 'NR>1 && toupper($1)==v {print $7}' "$OUT/fia_ba_sdi_percentiles.csv")
    printf 'group,sdimax\n%s,%s\n' "$var" "${sdi:-NA}" > "$OUT/sdimax_$v.csv"
    Rscript "$HERE/bakuzis_stress_v2.R" --input="$f" --source=projector --out="$OUT/score_projector/$v" --sdimax="$OUT/sdimax_$v.csv" --obs_ba_p999="${p999:-NA}" --label="$v" > "$OUT/score_projector_$v.log" 2>&1
    grep -E '^Reineke axis|Error' "$OUT/score_projector_$v.log" | sed "s/^/  $v: /"
  done
  echo "projector rescore written under $OUT/score_projector; the site axis is NOT EVALUATED unless --site_col names a column in the metrics file."
fi
echo "done. Nothing was pushed; close the job with github-manager when the outputs are reviewed."
