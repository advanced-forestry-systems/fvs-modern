#!/bin/bash
#SBATCH --job-name=bak_finish
#SBATCH --account=PUOM0008
#SBATCH --time=00:40:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/finish_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027
export PYTHONPATH=$W:$W/deployment/fvs2py:${PYTHONPATH:-}
cd $J
echo "== seed gate $(date)"; Rscript seedgate_full.R
echo "== aggregate $(date)"; mkdir -p ens_full_agg
python3 $W/calibration/python/bakuzis_uncertainty_aggregate.py --input-dir $J/ens_full --output-dir $J/ens_full_agg 2>&1 | tail -5
ls -la ens_full_agg
echo "== score $(date)"
Rscript -e 'suppressMessages(library(data.table)); p<-fread("proj_rescore/fia_ba_sdi_percentiles.csv"); o<-p[variant%in%c("NE"),.(group=variant,ba_p99,ba_p999)]; o<-rbind(o,o[,.(group="ACD",ba_p99,ba_p999)]); fwrite(o,"ens_obs_ba.csv"); s<-p[variant=="NE",.(group="NE",sdimax=round(sdi_p999_ha))]; s<-rbind(s,s[,.(group="ACD",sdimax)]); fwrite(s,"ens_sdimax.csv"); print(o); print(s)'
mkdir -p score_ens_full
Rscript bakuzis_stress_v2.R --input=ens_full_agg/bakuzis_uncertainty_summary_long.csv --source=uncert --out=score_ens_full --sdimax=ens_sdimax.csv --obs_ba_file=ens_obs_ba.csv --label=ens_full_20261028 > score_ens_full/run_report.txt 2>&1
grep -E '^Reineke axis|Error|Wrote' score_ens_full/run_report.txt
cat ens_full_agg/bakuzis_laws_compliance.csv 2>/dev/null
echo FINISH_DONE
