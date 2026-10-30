#!/bin/bash
## runoff_b_project.sh VAR ARM MODE BIN OUTROOT [VER]
## Contestant B: fvs-conus projector, VER = v2 (conus_eq_projector_v2.R, the Paper 2 Table 6 / proj_rescore2
## version) or v6 (conus_eq_projector_v6.R, v4 + HT-DBH monotonicity guard), on the catalog stands written by
## runoff_engine_runs.py --b-inputs. Ingrowth off (zero lookup) for both. Outputs under OUTROOT/VER/.
##   long: 20 cycles x 5 yr.  near: one cycle of k yr, one projector call per distinct k (K_CYCLE).
##   v2 takes --ncycles/--cyclelen; v6 takes --horizon/--cyclelen and needs --treelists to write tree lists.
## ARM: b1_native (mode free, logit mortality), b2_native (dependent, logit), b1_gompit, b2_gompit.
set -eu
VAR=$1; ARM=$2; MODE=$3; BIN=$4; OUT0=$5; VER=${6:-v2}; OUT=$OUT0/$VER
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj; P=$D/conus_eq_projector_${VER}.R
ZERO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/b_inputs/ingrowth_zero.rds
[ -f "$P" ] || { echo "MISSING $P"; exit 4; }
case $ARM in
  b1_native) MARGS="--mode=free";;
  b2_native) MARGS="--mode=dependent";;
  b1_gompit) MARGS="--mode=free --mort=gompit --cch_module=$D/cch_module.R";;
  b2_gompit) MARGS="--mode=dependent --mort=gompit --cch_module=$D/cch_module.R";;
  *) echo "bad arm $ARM"; exit 2;;
esac
step() { if [ "$VER" = v2 ]; then echo "--ncycles=$1 --cyclelen=$2"; else echo "--horizon=$(( $1 * $2 )) --cyclelen=$2 --treelists"; fi; }
v=$(echo $VAR | tr A-Z a-z); SI=$BIN/standinit_${VAR}.csv; TI=$BIN/treeinit_${v}
[ -f "$SI" ] || { echo "MISSING $SI"; exit 3; }
if [ "$MODE" = long ]; then
  O=$OUT/long/$ARM/$v; mkdir -p $O; cp -n $D/b1_re_means.rds $O/ 2>/dev/null || true
  Rscript $P --variant=$VAR $MARGS --nstands=0 $(step 20 5) --standinit_dir=$BIN --treeinit_dir=$TI --ingrowth=$ZERO --outdir=$O
else
  KS=$(Rscript -e "x<-read.csv('$SI'); k<-x\$K_CYCLE; cat(sort(unique(k[is.finite(k)])))")
  echo "distinct k: $KS"
  for k in $KS; do
    S=$OUT/near/$ARM/$v/k$k; mkdir -p $S/si; cp -n $D/b1_re_means.rds $S/ 2>/dev/null || true
    Rscript -e "x<-read.csv('$SI', colClasses=c(STAND_CN='character',STAND_ID='character')); write.csv(x[is.finite(x\$K_CYCLE) & x\$K_CYCLE==$k,], '$S/si/standinit_${VAR}.csv', row.names=FALSE)"
    Rscript $P --variant=$VAR $MARGS --nstands=0 $(step 1 $k) --standinit_dir=$S/si --treeinit_dir=$TI --ingrowth=$ZERO --outdir=$S
  done
fi
Rscript /fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/scripts/runoff_b_post.R $OUT0 $VAR $ARM $MODE $BIN $VER
echo "B_TASK_DONE $VER $VAR $ARM $MODE"
