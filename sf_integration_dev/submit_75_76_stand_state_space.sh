#!/bin/bash
# =============================================================================
# submit_75_76_stand_state_space.sh -- launch the García state-space transition
# fits on Cardinal: 75 top-height (H2|H1 GADA transition) and 76 stem-density
# (N(t) transition). Symmetric to submit_72_stand_bagrowth.sh. brms, 4 chains,
# seeded 20260702. Prints both job IDs; does NOT wait.
# Usage:  bash submit_75_76_stand_state_space.sh
# =============================================================================
set -euo pipefail
SLURM_ACCOUNT="${SLURM_ACCOUNT:-PUOM0008}"
WT="${WT:-/fs/scratch/PUOM0008/crsfaaron/wt-gompit}"
PAIRS="${PAIRS:-/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds}"
OUT_TH="${OUT_TH:-/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/stand_topheight}"
OUT_ST="${OUT_ST:-/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/stand_stems}"
mkdir -p "${OUT_TH}/logs" "${OUT_ST}/logs"

TH_JOB=$(sbatch --parsable <<SBATCH
#!/bin/bash
#SBATCH --job-name=stand_topheight
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=${OUT_TH}/logs/stand_topheight_%j.out
#SBATCH --error=${OUT_TH}/logs/stand_topheight_%j.err
set -euo pipefail
module load gcc/12.3.0 R/4.4.0
cd ${WT}/sf_integration_dev
Rscript 75_fit_stand_topheight.R \
  --pairs=${PAIRS} --out_dir=${OUT_TH} --n_sub=120000 --top_n=100
SBATCH
)
echo "SUBMITTED stand_topheight job: ${TH_JOB}"
echo "${TH_JOB}" > "${OUT_TH}/.last_job_id"

ST_JOB=$(sbatch --parsable <<SBATCH
#!/bin/bash
#SBATCH --job-name=stand_stems
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=${OUT_ST}/logs/stand_stems_%j.out
#SBATCH --error=${OUT_ST}/logs/stand_stems_%j.err
set -euo pipefail
module load gcc/12.3.0 R/4.4.0
cd ${WT}/sf_integration_dev
Rscript 76_fit_stand_stems.R \
  --pairs=${PAIRS} --out_dir=${OUT_ST} --n_sub=120000 --top_n=100 --offset_center=3.9
SBATCH
)
echo "SUBMITTED stand_stems job: ${ST_JOB}"
echo "${ST_JOB}" > "${OUT_ST}/.last_job_id"

echo "TOPHEIGHT_JOB=${TH_JOB}"
echo "STEMS_JOB=${ST_JOB}"
