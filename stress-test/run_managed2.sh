#!/bin/bash
SD=/fs/scratch/PUOM0008/crsfaaron/fvs_stress
module load python/3.12 2>/dev/null
for c in default calibrated gompit; do
  echo "### $c"
  python3 $SD/fvs_managed_v2.py --campaign $SD/out_fvs_v2 --config $c \
    --plantation $SD/plt_plantation.csv --rates $SD/state_harvest_rates.csv \
    --start 2030 --k 1.9 --window 20 --out $SD/managed2_$c
done
echo ALL_DONE
