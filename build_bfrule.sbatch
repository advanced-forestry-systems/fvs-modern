#!/bin/bash
#SBATCH --job-name=build_bfrule
#SBATCH --account=puom0008
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/logs/build_bfrule_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/logs/build_bfrule_%j.err

set -euo pipefail
cd /users/PUOM0008/crsfaaron/fvs-modern
module load gcc/12.3.0
bash deployment/scripts/build_fvs_executables.sh . /users/PUOM0008/crsfaaron/fvs-modern/bin-vol ne
echo "BUILD_EXIT=$?"
ls -la /users/PUOM0008/crsfaaron/fvs-modern/bin-vol/
nm /users/PUOM0008/crsfaaron/fvs-modern/bin-vol/FVSne | grep -i gregloadhg_ || echo "GREGLOADHG_SYMBOL_MISSING"
