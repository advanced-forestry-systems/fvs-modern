#!/bin/bash
module load gcc/12.3.0 python/3.12 >/dev/null 2>&1
cd ~/fvs-modern
export FVS_PROJECT_ROOT=$HOME/fvs-modern FVS_LIB_DIR=$HOME/fvs-modern/lib
SD=/fs/scratch/PUOM0008/crsfaaron/fvs_stress; TI=/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h
cp config/calibrated/ws.json /tmp/ws_orig.json
for V in 1.0 0.6 0.4; do
python3 - "$V" <<PY
import json,sys
val=float(sys.argv[1])
d=json.load(open("/tmp/ws_orig.json"))
d["calibration_multipliers"]["dds_multiplier"]=[val]*len(d["calibration_multipliers"]["dds_multiplier"])
json.dump(d,open("config/calibrated/ws.json","w"))
PY
python3 $SD/test_sdimax.py --variant ws --state CA --standinit-dir $SD/standinit_by_variant --treeinit-dir $TI --n-plots 40 --num-cycles 16 --tag "dds=$V"
done
cp /tmp/ws_orig.json config/calibrated/ws.json; echo restored
