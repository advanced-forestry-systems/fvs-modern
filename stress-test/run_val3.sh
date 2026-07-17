module load gcc/12.3.0 python/3.12 >/dev/null 2>&1
export FVS_PROJECT_ROOT=$HOME/fvs-modern FVS_LIB_DIR=$HOME/fvs-modern/lib FVS_FIA_DATA_DIR=/fs/scratch/PUOM0008/crsfaaron/FIA
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress
python3 run_gompit_projection.py --variant ne --n 6 \
  --coeff /fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv \
  --standinit-dir standinit_by_variant --treeinit-dir /fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit \
  --output out_gompit/ne_iter.csv --num-cycles 20 --cycle-length 5 --regen off \
  --arms default,calibrated,iterative > out_gompit/ne_iter.log 2>&1
echo "EXIT=$?" >> out_gompit/ne_iter.log
