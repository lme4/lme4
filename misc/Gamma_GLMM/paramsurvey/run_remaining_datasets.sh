#!/usr/bin/env bash
# Runs the R-side (glmmTMB/jointphi/lme4current) and Julia-side
# (MixedModels.jl pa/dispersion-again, via the warm DaemonMode server)
# fits for epil2_complex, report4bb, and schizophrenia -- the three
# paramsurvey datasets not yet covered (epil2_simple already done).
#
# Datasets are processed one at a time; within each dataset the three R
# methods run concurrently via mclapply (MC_CORES=10 each, 30 cores
# total) so we never oversubscribe past the requested 24-core budget.
# The Julia fit runs afterward through the already-warm daemon (fast,
# single client at a time -- not the bottleneck).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

MC_CORES=10
DATASETS=(report4bb schizophrenia epil2_complex)

for ex in "${DATASETS[@]}"; do
  echo "############################################"
  echo "### $ex : R-side fits (glmmTMB/jointphi/lme4current), MC_CORES=$MC_CORES each"
  echo "############################################"
  Rscript 02_fit_glmmTMB.R    "$ex" $MC_CORES > "run_${ex}_glmmTMB.log"    2>&1 &
  Rscript 03_fit_jointphi.R   "$ex" $MC_CORES > "run_${ex}_jointphi.log"   2>&1 &
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
echo "### ALL DATASETS COMPLETE"
echo "############################################"
