#!/bin/bash
#SBATCH --job-name=fvs_compile_ls_cs_ca_so
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:30:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/compile_ls_cs_ca_so_%j.out

module load gcc/12.3.0

echo "=== FVS variant compile: ls cs ca so ==="
echo "Started: $(date)"
echo "Host: $(hostname)"

cd /users/PUOM0008/crsfaaron/fvs-modern

env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_executables.sh \
    src-converted \
    src-converted/bin \
    ls cs ca so

echo "=== Build complete: $(date) ==="
echo "=== Final bin contents ==="
ls -lh src-converted/bin/FVS* | grep -v sourceList | grep -v "_sourceList"

