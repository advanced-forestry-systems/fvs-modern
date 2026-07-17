#!/bin/bash
module load gcc/12.3.0 python/3.12 >/dev/null 2>&1
cd ~/fvs-modern
export FVS_PROJECT_ROOT=$HOME/fvs-modern FVS_LIB_DIR=$HOME/fvs-modern/lib
SD=/fs/scratch/PUOM0008/crsfaaron/fvs_stress
TI=/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h
cp config/calibrated/ws.json /tmp/ws_orig.json
echo "RUN CURRENT (redwood NA->1052 default):"
python3 $SD/test_sdimax.py --variant ws --state CA --standinit-dir $SD/standinit_by_variant --treeinit-dir $TI --n-plots 40 --num-cycles 16 --tag CURRENT
cp config/calibrated_sdifix/ws.json config/calibrated/ws.json
echo "RUN SDIFIX (redwood NA->394 median):"
python3 $SD/test_sdimax.py --variant ws --state CA --standinit-dir $SD/standinit_by_variant --treeinit-dir $TI --n-plots 40 --num-cycles 16 --tag SDIFIX
cp /tmp/ws_orig.json config/calibrated/ws.json
echo "restored ws.json"
