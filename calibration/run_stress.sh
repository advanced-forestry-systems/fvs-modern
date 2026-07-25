#!/bin/bash
module purge
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
Rscript --vanilla stress_test_components.R
