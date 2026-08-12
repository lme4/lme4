#!/usr/bin/env bash
# SAS-side counterpart to run_full_survey.sh: export CSVs for PROC
# GLIMMIX, run RSPL + Laplace on all four real datasets, ingest results
# back into R, fold into the combined analysis. UNTESTED end-to-end --
# no local SAS available while writing this; see sas/fitlib.sas's header
# and the README's "SAS side" section for what's most likely to need
# fixing on first real run.
#
# Point SAS_EXE at your actual SAS batch executable/wrapper before
# running -- the invocation (flags, even the executable name) varies by
# install (SAS 9.4 Foundation vs. Viya vs. OnDemand); the default below
# (`sas -sysin ... -log ... -print ...`) is the standard Unix SAS 9.4
# batch-mode form, not a verified path on your system:
#   SAS_EXE=/path/to/sas ./run_sas_survey.sh
#
# SAS's process exit code is NOT a reliable success signal (it can come
# back non-zero even for a clean run with just warnings/notes) -- the
# SAS invocation below is deliberately exempted from set -e (`|| true`)
# and this script greps each .log for "^ERROR" instead of trusting $?.
# Every other step (the R export/ingest/analysis calls) keeps the normal
# fail-fast behaviour the rest of this project's orchestration scripts use.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

SAS_EXE="${SAS_EXE:-sas}"
DATASETS=(epil2_simple report4bb schizophrenia epil2_complex)

echo "############################################"
echo "### SAS PROC GLIMMIX survey (RSPL + Laplace), SAS_EXE=$SAS_EXE"
echo "############################################"

for ex in "${DATASETS[@]}"; do
  echo "############################################"
  echo "### $ex : export CSV for SAS"
  echo "############################################"
  Rscript 12_export_csv_for_sas.R "$ex" > "run_${ex}_export_sas.log" 2>&1
  tail -n 3 "run_${ex}_export_sas.log"

  echo "############################################"
  echo "### $ex : PROC GLIMMIX (RSPL + Laplace), via fit_${ex}.sas"
  echo "############################################"
  ( cd sas && "$SAS_EXE" -sysin "fit_${ex}.sas" -log "fit_${ex}.log" -print "fit_${ex}.lst" ) || true
  if grep -q "^ERROR" "sas/fit_${ex}.log"; then
    echo "!!! ERROR lines found in sas/fit_${ex}.log -- inspect before trusting sas/results_${ex}_*.csv !!!"
    grep "^ERROR" "sas/fit_${ex}.log"
  fi
  tail -n 30 "sas/fit_${ex}.log"

  echo "############################################"
  echo "### $ex : ingest SAS results + combined analysis"
  echo "############################################"
  for method in rspl laplace; do
    if [ -f "sas/results_${ex}_${method}.csv" ]; then
      Rscript 13_ingest_sas_results.R "$ex" "$method" > "run_${ex}_ingest_sas_${method}.log" 2>&1
      tail -n 5 "run_${ex}_ingest_sas_${method}.log"
    else
      echo "!!! sas/results_${ex}_${method}.csv not found, skipping ingest -- check sas/fit_${ex}.log !!!"
    fi
  done

  Rscript 07_analysis.R "$ex" > "run_${ex}_analysis_sas.log" 2>&1
  tail -n 25 "run_${ex}_analysis_sas.log"

  echo "=== $ex COMPLETE ==="
done

echo "############################################"
echo "### ALL DATASETS COMPLETE (SAS arm)"
echo "############################################"
