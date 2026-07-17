#!/bin/bash
#SBATCH --job-name=greghook_build_ne
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# Build a TEST FVSne binary with the native Greg DG/HG hook into the worktree
# bin-greghook/ (NOT the production bin path). Gates on the four Greg symbols.
module load gcc/12.3.0
WT=/fs/scratch/PUOM0008/crsfaaron/fvs-modern-greg-hook
cd "$WT"; mkdir -p logs bin-greghook
bash deployment/scripts/build_fvs_executables.sh . ./bin-greghook ne 2>&1 | tail -40
echo "=== SYMBOL GATE ==="
nm "$WT/bin-greghook/FVSne" 2>/dev/null | grep -iE "gregloaddg_|gregdgv_|gregloadhg_|greghgv_|gregclim" | sort
