#!/bin/bash
#SBATCH --job-name=bk_gregclim
#SBATCH --account=PUOM0008
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/bk_%x_%j.log
## usage: sbatch --job-name=<tag> --export=ALL,TAG=<tag>,DD0=..,TD=..,PPT=..,DD18=.. bakuzis_gregclim.slurm
set -uo pipefail
source /etc/profile.d/lmod.sh 2>/dev/null; module purge; module load gcc/12.3.0 R/4.5.2
WT=/users/PUOM0008/crsfaaron/fvs-modern-wt-greg2026
H=$WT/validation/bakuzis_stress_v2
PY=/fs/scratch/PUOM0008/crsfaaron/wt-engine/fvsvenv/bin/python
LIB=/fs/scratch/PUOM0008/crsfaaron/lib_greg2026/FVSne.so
OUT=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/bakuzis_greg2026_20261109
export FVS_PROJECT_ROOT=$WT FVS_CONFIG_DIR=$WT/config FVS_GREG_CONFDIR=$WT/config PYTHONNOUSERSITE=1
export PYTHONPATH="$WT:$WT/deployment/fvs2py:${PYTHONPATH:-}"
export FVS_GOMPIT=1 FVS_GOMPIT_COEF=$WT/config/greg_mortality_coefficients.csv
export FVS_GREGDG=1 FVS_GREGDG_COEF=$WT/config/greg_dg_coefficients.csv
export FVS_GREGHG=1 FVS_GREGHG_COEF=$WT/config/greg_hg_coefficients.csv
export FVS_GREGCRW=1 FVS_GREGCRW_COEF=$WT/config/greg_crown_change_coefficients.csv
export FVS_GREG_DD0=$DD0 FVS_GREG_TD=$TD FVS_GREG_PPT_SM=$PPT FVS_GREG_DD18=$DD18
env | grep -E '^FVS_(GOMPIT|GREG)' | sort > $OUT/env_$TAG.txt
$PY $H/engine_bakuzis_v2.py --wt $WT --lib $LIB --tag $TAG --out $OUT/engine_$TAG > $OUT/engine_$TAG.log 2>&1
grep -E 'GATE SUMMARY|BAKUZIS_ENGINE_V2_DONE' $OUT/engine_$TAG.log
Rscript $H/bakuzis_stress_v2.R --input=$OUT/engine_$TAG/engine_bakuzis_100yr_$TAG.csv --source=engine --out=$OUT/score_$TAG --sdimax=$OUT/sdimax_ne_engine_tph.csv --label=$TAG > $OUT/score_$TAG.log 2>&1
grep -E '^Reineke axis|Wrote outputs|Error' $OUT/score_$TAG.log
grep -h "PC1=" $OUT/engine_$TAG/run_${TAG}_BK01/* 2>/dev/null | head -2
echo ARM_DONE
