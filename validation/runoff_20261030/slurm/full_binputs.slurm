#!/bin/bash
#SBATCH --job-name=ro_full_binp --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full/logs/binp_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
RO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff; F=$RO/full; export FIA_DATA_DIR=$HOME/fia_data
mkdir -p $RO/b_inputs; [ -f $RO/b_inputs/ingrowth_zero.rds ] || Rscript $RO/scripts/make_zero_ingrowth.R $RO/b_inputs/ingrowth_zero.rds
for V in acd bm ca ci cr cs ec em ie ls ne pn sn so tt ut wc ws; do for M in long near; do
  python3 $RO/scripts/runoff_engine_runs.py --variant $V --arm native --mode $M --configs none --b-inputs $F/b_inputs/$M --out $F/b_inputs/$M/gate
done; done
echo BINP_DONE
