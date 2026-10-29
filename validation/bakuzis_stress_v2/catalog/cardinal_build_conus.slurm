#!/bin/bash
#SBATCH --job-name=fvs_build --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/build_conus_%j.out
# Build every CONUS variant library in the clean worktree after the GOMPSURV DBH-argument fix.
# The builder shadows base/PRGPRM.f90 in place, so builds are serial in one tree.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main; cd $W
echo "start $(date)"
# GOMPSURV fix: pass the tree DBH (D) as the sixth argument in every variant mortality routine.
cd src-converted
for f in $(grep -rl 'CALL GOMPSURV' --include=morts.f90 .); do
  sed -i -E 's/CALL GOMPSURV\(ISPC,CRG,(CCHG|CCHT\(I\)),(FINT|1\.0),SURVG\)/CALL GOMPSURV(ISPC,CRG,\1,\2,SURVG,D)/' $f
done
echo "calls without D: $(grep -rn 'CALL GOMPSURV' --include=morts.f90 . | grep -vc 'SURVG,D)')"
cd $W
mkdir -p lib_prefix; cp -n lib/FVSne.so lib_prefix/FVSne_prefix.so 2>/dev/null; cp -n lib/FVSacd.so lib_prefix/FVSacd_prefix.so 2>/dev/null
for v in acd bm ca ci cr cs ec em ie kt ls nc ne oc op pn so sn tt ut wc ws; do
  bash deployment/scripts/build_fvs_libraries.sh src-converted lib $v > build_$v.log 2>&1; echo "$v rc=$? $(date) $(ls -la lib/FVS$v.so 2>/dev/null | awk '{print $5}')"
done
ls -la lib/*.so
echo BUILD_CONUS_DONE
