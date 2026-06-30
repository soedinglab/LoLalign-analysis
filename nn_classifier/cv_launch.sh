#!/bin/bash
# Submit v20 (v17e + qsep6 + qfrac6 = 39 feats) 5-fold CV, chained after the qsep6+qfrac6 parquet build.
# Usage: BUILD=<jobid> bash cv_launch_v20.sh   (BUILD optional; if set, train waits on it)
set -euo pipefail
NN=/cbscratch/lolalignBenchmark/nn_classifier; cd $NN/slurm
DEP_TRAIN=""
if [[ -n "${BUILD:-}" ]]; then DEP_TRAIN="--dependency=afterok:$BUILD"; fi
FIN_IDS=()
for k in 0 1 2 3 4; do
  TAG=cv${k}_v20
  T=$(sbatch --parsable $DEP_TRAIN --export=ALL,TAG=$TAG cv_train.sh)
  P=$(sbatch --parsable --dependency=afterok:$T --export=ALL,TAG=$TAG cv_predict_array.sh)
  F=$(sbatch --parsable --dependency=afterok:$P --export=ALL,TAG=$TAG cv_finalize_one.sh)
  FIN_IDS+=("$F")
  echo "$TAG  train=$T predict=$P finalize=$F"
done
DEP=$(IFS=:; echo "${FIN_IDS[*]}")
AGG=$(sbatch --parsable --dependency=afterok:$DEP cv_aggregate_v20.sh)
echo "AGG=$AGG  (after finalizes: $DEP)"
echo "$AGG" > /tmp/cv_agg_v20_jobid.txt
