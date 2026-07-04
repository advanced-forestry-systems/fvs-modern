#!/bin/bash
#SBATCH --job-name=build_ne_exe
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/build_ne_exe_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/build_ne_exe_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
mkdir -p /fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs
echo "=== BUILD NE EXE START $(date) (with corrected gompmort) ==="
bash deployment/scripts/build_fvs_executables.sh . ./bin-test ne 2>&1 | tail -30
echo "=== RESULT ==="
ls -la /users/PUOM0008/crsfaaron/fvs-modern/bin-test/FVSne 2>/dev/null
nm /users/PUOM0008/crsfaaron/fvs-modern/bin-test/FVSne 2>/dev/null | grep -i gompsurv | head
echo "=== DONE $(date) ==="
