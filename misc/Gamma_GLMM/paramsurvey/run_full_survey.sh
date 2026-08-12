#!/usr/bin/env bash
# Full paramsurvey run at a given B (default 500) for all four datasets:
# prep (resimulate) -> CSV export for Julia -> R-side fits
# (glmmTMB/jointphi/lme4current, concurrent via mclapply, MC_CORES=10
# each = 30 cores total) -> Julia fit (MixedModels.jl pa/dispersion-
# again, via the warm DaemonMode server) -> ingest -> combined analysis.
# Datasets processed sequentially (fastest first) to stay within the
# 30-core budget; total wall time calibrated from the B=50 run is
# roughly 1.5-1.6 hours for B=500.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

B="${1:-500}"
MC_CORES=10
DATASETS=(epil2_simple report4bb schizophrenia epil2_complex)

echo "############################################"
echo "### FULL SURVEY, B=$B, MC_CORES=$MC_CORES per method (30 cores total)"
echo "############################################"

for ex in "${DATASETS[@]}"; do
  echo "############################################"
  echo "### $ex : prep (resimulate B=$B)"
  echo "############################################"
  Rscript "01_prep_${ex}.R" "$B" > "run_${ex}_prep.log" 2>&1
  tail -n 5 "run_${ex}_prep.log"

  echo "############################################"
  echo "### $ex : export CSV for Julia"
  echo "############################################"
  Rscript 09_export_csv_for_julia.R "$ex" > "run_${ex}_export.log" 2>&1
  tail -n 3 "run_${ex}_export.log"

  echo "############################################"
  echo "### $ex : R-side fits (glmmTMB/jointphi/lme4current), MC_CORES=$MC_CORES each"
  echo "############################################"
  Rscript 02_fit_glmmTMB.R     "$ex" $MC_CORES > "run_${ex}_glmmTMB.log"     2>&1 &
  Rscript 03_fit_jointphi.R    "$ex" $MC_CORES > "run_${ex}_jointphi.log"    2>&1 &
  Rscript 04_fit_lme4current.R "$ex" $MC_CORES > "run_${ex}_lme4current.log" 2>&1 &
  wait
  echo "--- $ex R-side done ---"
  tail -n 5 "run_${ex}_glmmTMB.log" "run_${ex}_jointphi.log" "run_${ex}_lme4current.log"

  echo "############################################"
  echo "### $ex : Julia fit (MixedModels.jl pa/dispersion-again)"
  echo "############################################"
  ( cd julia && ./run_via_daemon.sh "fit_${ex}.jl" ) > "run_${ex}_julia.log" 2>&1
  tail -n 5 "run_${ex}_julia.log"

  echo "############################################"
  echo "### $ex : ingest Julia results + combined analysis"
  echo "############################################"
  Rscript 10_ingest_julia_results.R "$ex" > "run_${ex}_ingest.log" 2>&1
  tail -n 5 "run_${ex}_ingest.log"
  Rscript 07_analysis.R "$ex" > "run_${ex}_analysis.log" 2>&1
  cat "run_${ex}_analysis.log"

  echo "=== $ex COMPLETE ==="
done

echo "############################################"
echo "### ALL DATASETS COMPLETE (B=$B)"
echo "############################################"
