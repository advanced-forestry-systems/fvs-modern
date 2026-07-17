module load gcc/12.3.0 python/3.12 >/dev/null 2>&1
export FVS_PROJECT_ROOT=$HOME/fvs-modern FVS_LIB_DIR=$HOME/fvs-modern/lib FVS_FIA_DATA_DIR=/fs/scratch/PUOM0008/crsfaaron/FIA
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress
python3 dbg2.py > dbg2.log 2>&1; echo DONE >> dbg2.log
