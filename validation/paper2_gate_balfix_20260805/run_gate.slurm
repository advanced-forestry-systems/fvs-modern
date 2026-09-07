#!/bin/bash
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --job-name=balfix_gate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805/gate_job_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805/gate_job_%j.err
set -u
G=/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805
R=/fs/scratch/PUOM0008/crsfaaron/balfix_run_20260804
mkdir -p $G
export TMPDIR=/fs/scratch/PUOM0008/crsfaaron/tmp

ND=$(ls $R/out_full/DONE_* 2>/dev/null | wc -l)
NS=$(ls $R/out_full/scorecard_fin_armA_s*.csv 2>/dev/null | wc -l)
echo "DONE files = $ND   scorecard shards = $NS"
if [ "$ND" -ne 32 ] || [ "$NS" -ne 32 ]; then
  echo "SHARDS INCOMPLETE. ABORTING GATE."
  echo "FAILED" > $G/GATE_STATUS
  exit 1
fi

echo "===== SCORECARD GATE ====="
python3 -u $G/gate_balfix.py       > $G/gate_balfix.log 2>&1
echo "gate_balfix.py exit=$?"
echo "===== TREELIST SANITY ====="
python3 -u $G/sanity_treelist.py   > $G/sanity_treelist.log 2>&1
echo "sanity_treelist.py exit=$?"
echo "OK" > $G/GATE_STATUS
touch $G/GATE_DONE
echo "end=$(date)"
