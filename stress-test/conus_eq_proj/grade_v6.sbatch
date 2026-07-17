#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J grade_v6
#SBATCH -t 0:40:00
#SBATCH -N 1 -n 4 --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/grade_v6_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=bakuzis_singlecohort
for b in v6_b2 v6_b2_dom; do
  Rscript build_singlecohort_traj_v6.R $b 10 40
  for sfx in _sc_all _sc_ne _sc_sn _sc_pn; do
    csv=$OUT/${b}${sfx}.csv
    [ -f "$csv" ] && Rscript bakuzis_check.R "$csv" 0.40 > $OUT/${b}${sfx}_check.txt 2>&1
    [ -f "$csv" ] && Rscript bakuzis_matrix.R "$csv" $OUT/${b}${sfx}_matrix.png 0.40 >/dev/null 2>&1 || true
  done
done
Rscript diag_stand.R > diag_stand_v6.txt 2>&1
for rg in NE SN PN; do Rscript --vanilla score_v5.R out_v6_seed/conus_eq_${rg,,}_conus_b2_metrics.csv $rg v6; done > score_v6_out.txt 2>&1
echo ALLDONE
