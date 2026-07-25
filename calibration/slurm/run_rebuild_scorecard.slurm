#!/bin/bash
#SBATCH --job-name=rebuild_scorecard
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/rebuild_scorecard_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/rebuild_scorecard_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

module load gcc/12.3.0
module load python/3.11
set -uo pipefail
TS=$(date +%Y%m%d_%H%M%S)
ROOT=/users/PUOM0008/crsfaaron/fvs-modern
echo "=== REBUILD+SCORECARD START $(date) node $(hostname) ==="

# --- Step 1: fix DG CSV double-# header (strip all leading comment lines) ---
DGCSV=/fs/scratch/PUOM0008/crsfaaron/wt-ne-dg/config/greg_dg_coefficients.csv
cp "$DGCSV" "${DGCSV}.bak_${TS}"
grep -v '^#' "${DGCSV}.bak_${TS}" > "$DGCSV"
echo "[DG csv] comment lines now: $(grep -c '^#' "$DGCSV")  first line: $(head -1 "$DGCSV")  data rows: $(( $(wc -l < "$DGCSV") - 1 ))"

# --- Step 2: clean rebuild of FVSne (fresh temp objects; source list has HG) ---
BINOUT=$ROOT/bin-rebuilt_${TS}
mkdir -p "$BINOUT"
cd "$ROOT"
echo "[build] starting"
bash deployment/scripts/build_fvs_executables.sh . "$BINOUT" ne 2>&1 | tail -25
NEWBIN="$BINOUT/FVSne"
if [ ! -x "$NEWBIN" ]; then echo "GATE FAIL: no binary produced at $NEWBIN"; exit 1; fi
ls -lh "$NEWBIN"

# --- Step 3: GATE on Greg HG object being linked ---
echo "[gate] checking symbols"
HG_SYM=$(nm "$NEWBIN" 2>/dev/null | grep -c "gregloadhg_")
HG_STR=$(strings "$NEWBIN" | grep -c "GREGHG enabled")
DG_SYM=$(nm "$NEWBIN" 2>/dev/null | grep -c "gregloaddg_")
echo "[gate] gregloadhg_=$HG_SYM  'GREGHG enabled'str=$HG_STR  gregloaddg_=$DG_SYM"
if [ "$HG_SYM" -lt 1 ]; then echo "GATE FAIL: Greg HG object still not linked. Aborting before scorecard."; exit 2; fi
echo "[gate] PASS: HG object present"

# --- Step 4: install rebuilt binary into runner path (backup old) ---
RUNBIN=$ROOT/src-converted/bin/FVSne
cp "$RUNBIN" "${RUNBIN}.pre_rebuild_${TS}.bak"
cp "$NEWBIN" "$RUNBIN"
echo "[install] runner binary md5: $(md5sum "$RUNBIN" | cut -d' ' -f1)  (backup ${RUNBIN}.pre_rebuild_${TS}.bak)"

# --- Step 5: rerun scorecard into a fresh dir ---
OUT=/fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_rebuilt
mkdir -p "$OUT"
cd "$ROOT/calibration"
echo "=== SCORECARD START $(date) ==="
python3 run_5model_scorecard.py --variant ne --max-pairs 5000 --outdir "$OUT" --fia-dir /fs/scratch/PUOM0008/crsfaaron/FIA
echo "[scorecard] exit=$?"

# --- Step 6: quick separation check (organon vs conus_climate) ---
CSV="$OUT/scorecard_ne_5model.csv"
if [ -f "$CSV" ]; then
python3 - "$CSV" << 'PY'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1])))
def g(m,c): return [float(r[c]) for r in rows if r['model']==m and r[c] not in ('','NA')]
import statistics as st
for m in ['fvs_base','organon','conus_spdep','conus_climate']:
    ba=g(m,'BA_PRED'); print(f"{m:14} n={len(ba)} meanBA_PRED={st.mean(ba):.3f}" if ba else f"{m}: none")
o=[r['BA_PRED'] for r in rows if r['model']=='organon']
c=[r['BA_PRED'] for r in rows if r['model']=='conus_climate']
maxd=max(abs(float(a)-float(b)) for a,b in zip(o,c)) if o and c else -1
print(f"max |organon-conus_climate| BA_PRED = {maxd}  (>0 means they now DIFFER, defect fixed)")
PY
fi
echo "=== REBUILD+SCORECARD DONE $(date) ==="
