#!/bin/bash
#SBATCH --job-name=gompsize_build
#SBATCH --account=PUOM0008
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/build_gompsize_%j.log
set -euo pipefail
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
echo "START $(date)"
gfortran --version | head -1
bash deployment/scripts/build_fvs_executables.sh . /users/PUOM0008/crsfaaron/fvs-modern/bin-gompsize ne
echo "BUILD_RC=$?"
echo "=== GATE ==="
ls -la /users/PUOM0008/crsfaaron/fvs-modern/bin-gompsize/ || true
B=/users/PUOM0008/crsfaaron/fvs-modern/bin-gompsize/FVSne
if [ -f "$B" ]; then
  echo "BINARY_EXISTS yes"
  nm "$B" 2>/dev/null | grep -iE "gompsurv_|gompload_|gregloadhg_" || echo "no matching symbols"
else
  echo "BINARY_EXISTS no"
fi
echo "END $(date)"
