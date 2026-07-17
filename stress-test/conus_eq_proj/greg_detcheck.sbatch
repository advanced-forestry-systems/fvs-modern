#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:40:00
#SBATCH -N 1 -n 8
#SBATCH --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/greg_detcheck_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
TMP=$D/detcheck; mkdir -p $TMP/seq $TMP/par
echo "=== SEQUENTIAL: NE 300-stand slice (seed 7) ==="
Rscript --vanilla $D/conus_eq_projector_greg.R --variant=NE --nstands=300 --seed=7 --outdir=$TMP/seq
echo "=== PARALLEL (8 cores): NE 300-stand slice (seed 7) ==="
Rscript --vanilla $D/conus_eq_projector_greg_par.R --variant=NE --nstands=300 --seed=7 --cores=8 --outdir=$TMP/par
echo "=== DIFF metrics ==="
S=$TMP/seq/conus_eq_ne_greg_metrics.csv; P=$TMP/par/conus_eq_ne_greg_metrics.csv
# sort both (parallel may reorder within state if mc.preschedule chunks differ) then compare
sort $S > $TMP/s_sorted.csv; sort $P > $TMP/p_sorted.csv
if diff -q $TMP/s_sorted.csv $TMP/p_sorted.csv >/dev/null; then echo "METRICS_IDENTICAL (sorted)"; else echo "METRICS_DIFFER"; diff $TMP/s_sorted.csv $TMP/p_sorted.csv | head -20; fi
echo "rows seq=$(wc -l <$S) par=$(wc -l <$P)"
# also raw (unsorted) row-order check
if diff -q $S $P >/dev/null; then echo "METRICS_IDENTICAL_RAWORDER"; else echo "RAWORDER_DIFFERS (row order only; sorted match is what matters)"; fi
echo "=== DIFF fallback ==="
sort $TMP/seq/conus_eq_ne_greg_fallback.csv > $TMP/sf.csv; sort $TMP/par/conus_eq_ne_greg_fallback.csv > $TMP/pf.csv
diff -q $TMP/sf.csv $TMP/pf.csv && echo "FALLBACK_IDENTICAL" || echo "FALLBACK_DIFFER"
echo "DETCHECK_DONE"
