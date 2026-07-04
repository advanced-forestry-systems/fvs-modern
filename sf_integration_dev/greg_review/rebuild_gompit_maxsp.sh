#!/bin/bash
set -e
cd /fs/scratch/PUOM0008/crsfaaron/wt-engine
module load gcc/12.3.0 2>/dev/null || true
# Shadow base/PRGPRM.f90 with ne MAXSP=108 for a consistent COMMON layout
cp -n src-converted/base/PRGPRM.f90 src-converted/base/PRGPRM.f90.orig
cp src-converted/ne/common/PRGPRM.f90 src-converted/base/PRGPRM.f90
echo "[rebuild] base/PRGPRM MAXSP now:" $(grep 'MAXSP *=' src-converted/base/PRGPRM.f90 | head -1)
rm -rf lib_gompit_fixed && mkdir -p lib_gompit_fixed
bash deployment/scripts/build_fvs_libraries.sh src-converted lib_gompit_fixed ne
RC=$?
# Restore original base/PRGPRM.f90
cp src-converted/base/PRGPRM.f90.orig src-converted/base/PRGPRM.f90
echo "[rebuild] restored base/PRGPRM MAXSP:" $(grep 'MAXSP *=' src-converted/base/PRGPRM.f90 | head -1)
echo "[rebuild] exit=$RC lib:"; ls -la lib_gompit_fixed/FVSne.so 2>/dev/null
echo REBUILD_DONE
