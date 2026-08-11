#!/bin/bash
# =============================================================================
# submit_72_stand_bagrowth.sh -- launch 72_fit_stand_bagrowth.R on Cardinal.
# Symmetric to the stand-survival fit (71). brms, 4 chains, seeded 20260702.
# Usage:  bash submit_72_stand_bagrowth.sh
# =============================================================================
set -euo pipefail
SLURM_ACCOUNT="${SLURM_ACCOUNT:-PUOM0008}"
WT="${WT:-/fs/scratch/PUOM0008/crsfaaron/wt-gompit}"
OUT_DIR="${OUT_DIR:-/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/stand_bagrowth}"
PAIRS="${PAIRS:-/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds}"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

JOB_ID=$(sbatch --parsable <<SBATCH
#!/bin/bash
#SBATCH --job-name=stand_bagrowth
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=${LOG_DIR}/stand_bagrowth_%j.out
#SBATCH --error=${LOG_DIR}/stand_bagrowth_%j.err

set -euo pipefail
module load gcc/12.3.0 R/4.4.0
cd ${WT}/sf_integration_dev
Rscript 72_fit_stand_bagrowth.R \
  --pairs=${PAIRS} \
  --out_dir=${OUT_DIR} \
  --n_sub=120000 --tau_m=10 --tau_d=15
SBATCH
)
echo "SUBMITTED stand_bagrowth job: ${JOB_ID}"
echo "${JOB_ID}" > "${OUT_DIR}/.last_job_id"
