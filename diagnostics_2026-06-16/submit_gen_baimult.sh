#!/bin/bash
#SBATCH --job-name=gen_baimult
#SBATCH --time=05:00:00
#SBATCH --mem=14G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/gen_baimult_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/gen_baimult_%j.err
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
declare -A VS=( [acd]="ME,NH,VT" [sn]="AL,GA,SC,MS" [ls]="MI,MN,WI" [cs]="IL,IN,MO" [ie]="ID,MT,WA" [kt]="MT,ID" [ci]="ID" [em]="MT,ND,SD" [bm]="OR,WA" [cr]="CO,WY" [tt]="WY,ID" [ut]="UT,NV" [ca]="CA" [ws]="CA" [nc]="CA,OR" [so]="CA,OR" [ec]="OR,WA" [wc]="OR,WA" [oc]="OR,WA" [op]="WA" [pn]="OR,WA" )
for v in "${!VS[@]}"; do
  echo "=== baimult $v (${VS[$v]}) ==="
  VAR=$v STATES="${VS[$v]}" NSAMP=150 SEED=5 python3 calib_ne.py 2>&1 | grep -E "calibrated species|DONE_CALIB|Error|Traceback" | head -3
done
echo ALL_BAIMULT_DONE
