#!/bin/bash
#SBATCH --job-name=crw_stress_fix
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0:20:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/crw_stress_%j.out
module load gcc/12.3.0
module load R/4.4.0
Rscript /users/PUOM0008/crsfaaron/fvs-modern/calibration/stress_test_components.R 2>&1
echo "Done: $(date)"
ls -lh /users/PUOM0008/crsfaaron/fvs-modern/calibration/output/stress_test/

