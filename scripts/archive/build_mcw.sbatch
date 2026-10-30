#!/bin/bash
#SBATCH --job-name=fvs_mcw_build
#SBATCH --account=PUOM0008
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/build_mcw_%j.log
set -euo pipefail
cd /users/PUOM0008/crsfaaron/fvs-modern
module load gcc/12.3.0
echo "=== build start $(date) ==="
bash deployment/scripts/build_fvs_executables.sh . /users/PUOM0008/crsfaaron/fvs-modern/bin-mcw ne
echo "=== build end $(date) ==="
echo "=== gate ==="
ls -la /users/PUOM0008/crsfaaron/fvs-modern/bin-mcw/FVSne
nm /users/PUOM0008/crsfaaron/fvs-modern/bin-mcw/FVSne | grep gregloadhg_ || echo "GATE_FAIL_no_gregloadhg"
echo "=== DONE ==="
