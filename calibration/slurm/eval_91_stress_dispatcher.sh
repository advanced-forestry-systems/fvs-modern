#!/bin/bash
# =============================================================================
# Comprehensive site productivity stress test dispatcher
# =============================================================================
# Submits one Stan refit per (model, variant, site_var) combination. Stress
# tests answer the manuscript question: are the model conclusions robust to
# the choice of site productivity variable (cspi vs bgi vs climate_si)?
#
# Coverage matrix:
#   model      | B1 species-free | B2 species-aware | site variables
#   -----------|-----------------|------------------|-----------------
#   HG         | YES (binary)    | YES (banked Apr25)| cspi*, bgi, climate_si
#   DG Kuehne  | YES (May 13)    | YES (May 7)       | cspi*, bgi, climate_si
#   HT-DBH     | NO (architecture)| YES (Apr 29)     | cspi*, bgi, climate_si
#   CR         | NO              | YES (Apr 29)      | cspi*, bgi, climate_si
#   HCB        | NO              | YES (May 3)       | cspi*, bgi, climate_si
#   Mortality  | NO              | YES (Apr 29)      | cspi*, bgi, climate_si
#
# * = cspi already banked, no new fit needed
#
# Total new fits if all submitted: 8 (HG+DG_Kue B1+B2 × bgi+climate_si)
# Plus per-model B2 × bgi+climate_si for the other four: 8 more
# Grand total comprehensive: 16 fits at ~$80-120 each = ~$1,500-1,900
#
# Usage:
#   bash eval_91_stress_test_dispatcher.sh hg_dg            # priority subset
#   bash eval_91_stress_test_dispatcher.sh all_b2           # B2 stress for all 6 models
#   bash eval_91_stress_test_dispatcher.sh full             # everything
# =============================================================================

set -uo pipefail

PROJ_ROOT="/users/PUOM0008/crsfaaron/fvs-modern"
cd "$PROJ_ROOT"

submit_site_stress() {
  local MODEL=$1
  local VARIANT=$2
  local SITE_VAR=$3
  local SUBSAMPLE=${4:-200000}
  local WALLTIME=${5:-18:00:00}

  cat > /tmp/stress_${MODEL}_${VARIANT}_${SITE_VAR}.sh << EOF
#!/bin/bash
#SBATCH --job-name=stress_${MODEL}_${VARIANT}_${SITE_VAR}
#SBATCH --time=${WALLTIME}
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=${PROJ_ROOT}/logs/stress_${MODEL}_${VARIANT}_${SITE_VAR}_%j.out
#SBATCH --error=${PROJ_ROOT}/logs/stress_${MODEL}_${VARIANT}_${SITE_VAR}_%j.err

module load gcc/12.3.0 R/4.4.0
cd ${PROJ_ROOT}

Rscript --vanilla calibration/R/eval/90_site_productivity_stress_test.R \\
  --model=${MODEL} \\
  --variant=${VARIANT} \\
  --site_var=${SITE_VAR} \\
  --subsample=${SUBSAMPLE} \\
  --outdir=calibration/output/conus/${MODEL}/site_stress
EOF
  echo ">> sbatch stress_${MODEL}_${VARIANT}_${SITE_VAR}"
  sbatch /tmp/stress_${MODEL}_${VARIANT}_${SITE_VAR}.sh
}

dispatch_hg_dg() {
  # Priority: HG and DG_Kue with both variants × non-cspi site variables
  # (cspi versions already banked)
  for SITE_VAR in bgi climate_si; do
    for spec in "dg_kue b1 200000" "dg_kue b2 200000" \
                "hg b1 500000" "hg b2 500000"; do
      read MODEL VARIANT SUBSAMPLE <<< $spec
      submit_site_stress $MODEL $VARIANT $SITE_VAR $SUBSAMPLE
    done
  done
}

dispatch_all_b2() {
  # B2 stress test for the four other base models × non-cspi site variables
  for SITE_VAR in bgi climate_si; do
    for spec in "htdbh b2 500000 24:00:00" \
                "cr b2 500000 24:00:00" \
                "hcb b2 500000 24:00:00" \
                "mort b2 500000 24:00:00"; do
      read MODEL VARIANT SUBSAMPLE WALL <<< $spec
      submit_site_stress $MODEL $VARIANT $SITE_VAR $SUBSAMPLE $WALL
    done
  done
}

case "${1:-hg_dg}" in
  hg_dg)
    dispatch_hg_dg
    ;;
  all_b2)
    dispatch_all_b2
    ;;
  full)
    dispatch_hg_dg
    dispatch_all_b2
    ;;
  *)
    echo "Usage: $0 {hg_dg|all_b2|full}"
    exit 1
    ;;
esac

echo ""
echo "After fits land, run the comparison:"
echo "  for M in hg dg_kue htdbh cr hcb mort; do"
echo "    Rscript calibration/R/eval/90b_site_stress_compare.R --model=\$M"
echo "  done"
