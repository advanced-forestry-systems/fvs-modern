#!/bin/bash
# Run-off full run (A near, C long+near, B v2+v6 long+near), then scoring and decision.
cd /fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full
J1=$(sbatch --parsable full_engines.slurm)
J2=$(sbatch --parsable full_binputs.slurm)
J3=$(sbatch --parsable --dependency=afterok:$J2 full_b.slurm)
J4=$(sbatch --parsable --dependency=afterany:$J1:$J3 full_score.slurm)
echo "$(date '+%F %T') engines=$J1 binputs=$J2 b=$J3 score=$J4" | tee -a submitted.txt
