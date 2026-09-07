#!/bin/bash
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --job-name=bal3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805/bal3_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805/bal3_%j.err
set -u
G=/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805
export TMPDIR=/fs/scratch/PUOM0008/crsfaaron/tmp
cd $G
python3 -u sanity_bal.py > $G/sanity_bal.log 2>&1
echo "exit=$?"
touch $G/BAL3_DONE
