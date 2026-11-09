#!/bin/bash
#SBATCH --job-name=bakuzis_greg2026
#SBATCH --account=PUOM0008
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/bakuzis_greg2026_%j.log
set -uo pipefail
source /etc/profile.d/lmod.sh 2>/dev/null; module purge; module load gcc/12.3.0 R/4.5.2
WT=/users/PUOM0008/crsfaaron/fvs-modern-wt-greg2026
H=$WT/validation/bakuzis_stress_v2
PY=/fs/scratch/PUOM0008/crsfaaron/wt-engine/fvsvenv/bin/python
LIB=/fs/scratch/PUOM0008/crsfaaron/lib_greg2026/FVSne.so
OUT=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/bakuzis_greg2026_20261109
IMP=/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact
mkdir -p $OUT
export FVS_PROJECT_ROOT=$WT FVS_CONFIG_DIR=$WT/config FVS_GREG_CONFDIR=$WT/config PYTHONNOUSERSITE=1
export PYTHONPATH="$WT:$WT/deployment/fvs2py:${PYTHONPATH:-}"
printf 'group,sdimax\nSpruce-Fir,1250\nNorthern-Hardwood,1384\nPine,1307\nOak-Pine,1384\n' > $OUT/sdimax_ne_engine_tph.csv
run_arm () {   # tag, then env assignments
  tag=$1; shift
  ( for kv in "$@"; do export "$kv"; done
    env | grep -E '^FVS_(GOMPIT|GREG)' | sort > $OUT/env_$tag.txt
    $PY $H/engine_bakuzis_v2.py --wt $WT --lib $LIB --tag $tag --out $OUT/engine_$tag > $OUT/engine_$tag.log 2>&1 )
  grep -E 'GATE SUMMARY|BAKUZIS_ENGINE_V2_DONE' $OUT/engine_$tag.log || echo "ARM $tag did not finish"
  f=$OUT/engine_$tag/engine_bakuzis_100yr_$tag.csv
  [ -f $f ] && Rscript $H/bakuzis_stress_v2.R --input=$f --source=engine --out=$OUT/score_$tag --sdimax=$OUT/sdimax_ne_engine_tph.csv --label=$tag > $OUT/score_$tag.log 2>&1
  grep -E '^Reineke axis|Wrote outputs|Error' $OUT/score_$tag.log
}
CLIM="FVS_GREG_DD0=1600 FVS_GREG_TD=28 FVS_GREG_PPT_SM=230 FVS_GREG_DD18=3500"
run_arm native
run_arm mortprev FVS_GOMPIT=1 FVS_GOMPIT_COEF=$IMP/greg_mortality_coefficients_main_fa29b44.csv
run_arm mortnew  FVS_GOMPIT=1 FVS_GOMPIT_COEF=$WT/config/greg_mortality_coefficients.csv
run_arm gregfull FVS_GOMPIT=1 FVS_GOMPIT_COEF=$WT/config/greg_mortality_coefficients.csv \
  FVS_GREGDG=1 FVS_GREGDG_COEF=$WT/config/greg_dg_coefficients.csv \
  FVS_GREGHG=1 FVS_GREGHG_COEF=$WT/config/greg_hg_coefficients.csv \
  FVS_GREGCRW=1 FVS_GREGCRW_COEF=$WT/config/greg_crown_change_coefficients.csv $CLIM
grep -h "GREGHG enabled\|GREGDG\|GOMP" $OUT/engine_gregfull/*.out 2>/dev/null | sort | uniq -c | head
echo ALL_DONE
