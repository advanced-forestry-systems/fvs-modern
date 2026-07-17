#!/bin/bash
SD=/fs/scratch/PUOM0008/crsfaaron/fvs_stress
module load python/3.12 2>/dev/null
python3 $SD/fvs_perseus_aggregate.py --in-dir $SD/out_fvs_v2 --out-dir $SD/perseus_series_default_v2 --config default --start 2025 --engine fvs_default
python3 $SD/fvs_perseus_aggregate.py --in-dir $SD/out_fvs_v2 --out-dir $SD/perseus_series_calibrated_v2 --config calibrated --start 2025 --engine fvs_calibrated
python3 $SD/fvs_perseus_aggregate.py --in-dir $SD/out_gompit_v2 --out-dir $SD/perseus_series_gompit_v2 --config gompit --start 2025 --engine fvs_gompit
echo ALL_AGG_DONE
