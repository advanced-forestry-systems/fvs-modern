#!/bin/bash
#SBATCH --job-name=bk_treelist
#SBATCH --account=PUOM0008
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/ooda_20261110/treelist_b0_%j.log
source /etc/profile.d/lmod.sh 2>/dev/null; module purge; module load gcc/12.3.0
WT=/users/PUOM0008/crsfaaron/fvs-modern-wt-greg2026
PY=/fs/scratch/PUOM0008/crsfaaron/wt-engine/fvsvenv/bin/python
LIB=/fs/scratch/PUOM0008/crsfaaron/lib_greg2026/FVSne.so
OUT=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/bakuzis_greg2026_20261109/treelist_20261110
D=/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/ooda_20261110
C=/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/climate/maine_greg_pc_climate_summary_1991-2020.csv
export FVS_PROJECT_ROOT=$WT FVS_CONFIG_DIR=$WT/config FVS_GREG_CONFDIR=$WT/config PYTHONNOUSERSITE=1 PYTHONPATH="$WT:$WT/deployment/fvs2py"
ALL="Northern-Hardwood,Oak-Pine,Spruce-Fir,Pine"; S="Low,Medium,High"
: native done
for g in Northern-Hardwood Oak-Pine Spruce-Fir Pine; do
  IFS=, read gg n s dd0 td ppt dd18 p1 p2 <<<"$(grep "^$g,.*,median," $C)"
  ( export FVS_GOMPIT=1 FVS_GOMPIT_COEF=$WT/config/greg_mortality_coefficients.csv FVS_GREGDG=1 FVS_GREGDG_COEF=$WT/config/greg_dg_coefficients.csv \
      FVS_GREGHG=1 FVS_GREGHG_COEF=$D/greg_hg_coefficients_B0ne999.csv FVS_GREGCRW=1 FVS_GREGCRW_COEF=$WT/config/greg_crown_change_coefficients.csv \
      FVS_GREG_DD0=$dd0 FVS_GREG_TD=$td FVS_GREG_PPT_SM=$ppt FVS_GREG_DD18=$dd18
    $PY $D/treelist_runs.py --wt $WT --lib $LIB --tag gregbne --out $OUT --groups $g --sites $S --dens $S )
done
echo ALL_DONE
