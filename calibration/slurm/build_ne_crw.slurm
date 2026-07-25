#!/bin/bash
#SBATCH --job-name=build_ne_crw
#SBATCH --account=PUOM0008
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/build_ne_crw_%j.log
module load gcc/12.3.0
echo "=== BUILD FVSne (GREGCRWV env-var fix) $(date) ==="
cd /users/PUOM0008/crsfaaron/fvs-modern
BTMP=/tmp/fvs-crw-$$
mkdir -p "$BTMP"
bash deployment/scripts/build_fvs_executables.sh . "$BTMP" ne 2>&1 | tail -30
if [ ! -f "$BTMP/FVSne" ]; then echo "BUILD FAILED"; exit 1; fi
cp "$BTMP/FVSne" src-converted/bin/FVSne
rm -rf "$BTMP"
ls -la src-converted/bin/FVSne
echo "=== DONE $(date) ==="
