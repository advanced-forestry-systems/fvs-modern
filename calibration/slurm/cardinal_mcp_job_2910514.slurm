#!/bin/bash
#SBATCH --job-name=fvs_rebuild_all
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/rebuild_all_%j.out

module load gcc/12.3.0

cd /users/PUOM0008/crsfaaron/fvs-modern/src-converted

# Touch patched file to force recompile of ne; touch others to be consistent
touch ne/dgf.f90 sn/dgf.f90 wc/dgf.f90 em/dgf.f90 pn/dgf.f90

# Build each variant
for var in ne sn wc em pn; do
  echo "=== Building FVS${var} ==="
  bash deployment/scripts/build_fvs_executables.sh $var 2>&1 | tail -5
  echo "FVS${var} exit: $?"
done

echo ""
echo "=== Final binaries ==="
ls -lh bin/FVSne bin/FVSsn bin/FVSwc bin/FVSem bin/FVSpn

