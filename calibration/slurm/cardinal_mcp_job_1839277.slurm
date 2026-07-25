#!/bin/bash
#SBATCH --job-name=rebuild_fvsne
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/rebuild_ne_%j.out

module load gcc/12.3.0

cd /users/PUOM0008/crsfaaron/fvs-modern

echo "=== Building FVSne ==="
echo "Started: $(date)"
echo "Host: $(hostname)"

bash /users/PUOM0008/crsfaaron/fvs-modern/deployment/scripts/build_fvs_executables.sh \
    /users/PUOM0008/crsfaaron/fvs-modern/src-converted \
    /users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin \
    ne

echo "Build exit: $?"
echo "=== Done: $(date) ==="
ls -lh /users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin/FVSne

