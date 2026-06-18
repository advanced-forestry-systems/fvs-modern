#!/bin/bash
#SBATCH --job-name=gdb_inproc
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/gdb_inproc_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/gdb_inproc_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export LIB=/users/PUOM0008/crsfaaron/fvs-modern/lib-test/FVSne.so
gdb -batch -ex 'set pagination off' -ex run -ex 'bt' -ex 'info registers rip' -ex 'quit' --args python3 test_inproc_diag.py 2>&1 | tail -45
