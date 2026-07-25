#!/bin/bash
#SBATCH --job-name=dg_pc_deploy
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/dg_pc_deploy_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/dg_pc_deploy_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

module load gcc/12.3.0

echo "=== DG PC deploy started: $(date) ==="

cd /users/PUOM0008/crsfaaron/
python3 /users/PUOM0008/crsfaaron/cch-fix-bundle/deploy_dg_pc_update.py 2>&1 | tee /tmp/pc_deploy.log

DEPLOY_RC=${PIPESTATUS[0]}
echo "=== deploy_dg_pc_update.py exited with RC=${DEPLOY_RC} ==="

if [ ${DEPLOY_RC} -eq 0 ]; then
  echo "=== Verifying binaries ==="
  for VARIANT in ne sn wc em pn; do
    WT="/fs/scratch/PUOM0008/crsfaaron/wt-${VARIANT}-dg"
    EXE="${WT}/src-converted/bin/FVS${VARIANT}"
    if [ -f "${EXE}" ]; then
      MTIME=$(stat -c '%y' "${EXE}")
      SYMS=$(nm "${EXE}" 2>/dev/null | grep -icE "gregdgv_|gregloaddg_" || echo 0)
      echo "  FVS${VARIANT}: exists mtime=${MTIME} greg_syms=${SYMS}"
    else
      echo "  FVS${VARIANT}: binary NOT FOUND at ${EXE}"
    fi
  done
fi

echo "=== Deploy job complete: $(date) ==="
exit ${DEPLOY_RC}
