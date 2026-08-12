#!/usr/bin/env bash
# Fill in what's needed to regenerate the original six-method summary
# plot (08_summary_plots.R) at B=500: build the isolated lme4 2.0-6
# library (one-time, not committed -- see 05_fit_lme4old.R), prep
# epil2_phigt1 (never run at B=500 this session, has no results at all
# yet), fit glmmTMB/jointphi/lme4current for epil2_phigt1 specifically,
# then fit pirlsdigamma/pirlsmoment/lme4old for all five datasets (only
# glmmTMB/jointphi/lme4current/juliaMixedModels exist so far for the
# other four, from the Julia-comparison run -- reuse their existing
# B=500 simdata as-is, don't resimulate).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

B=500
MC_CORES=10

echo "############################################"
echo "### build isolated lme4 2.0-6 library"
echo "############################################"
if [ ! -d lme4_206_lib ]; then
  Rscript -e '
    tmplib <- here::here("misc/Gamma_GLMM/paramsurvey/lme4_206_lib")
    dir.create(tmplib, showWarnings = FALSE)
    withr::with_libpaths(tmplib, remotes::install_version("lme4", "2.0-6", upgrade = "never"))
  ' > run_build_lme4_206_lib.log 2>&1
  tail -n 20 run_build_lme4_206_lib.log
else
  echo "lme4_206_lib already exists, skipping"
fi

echo "############################################"
echo "### epil2_phigt1 : prep (B=$B) -- never run at B=500 before"
echo "############################################"
Rscript 01_prep_epil2_phigt1.R "$B" > run_epil2_phigt1_prep.log 2>&1
tail -n 5 run_epil2_phigt1_prep.log

echo "############################################"
echo "### epil2_phigt1 : glmmTMB/jointphi/lme4current, MC_CORES=$MC_CORES each"
echo "############################################"
Rscript 02_fit_glmmTMB.R     epil2_phigt1 $MC_CORES > run_epil2_phigt1_glmmTMB.log     2>&1 &
Rscript 03_fit_jointphi.R    epil2_phigt1 $MC_CORES > run_epil2_phigt1_jointphi.log    2>&1 &
Rscript 04_fit_lme4current.R epil2_phigt1 $MC_CORES > run_epil2_phigt1_lme4current.log 2>&1 &
wait
echo "--- epil2_phigt1 core three done ---"
tail -n 5 run_epil2_phigt1_glmmTMB.log run_epil2_phigt1_jointphi.log run_epil2_phigt1_lme4current.log

DATASETS=(epil2_simple epil2_phigt1 epil2_complex report4bb schizophrenia)

for ex in "${DATASETS[@]}"; do
  echo "############################################"
  echo "### $ex : pirlsdigamma/pirlsmoment/lme4old, MC_CORES=$MC_CORES each"
  echo "############################################"
  Rscript 06_fit_pirls_phi.R "$ex" digamma $MC_CORES > "run_${ex}_pirlsdigamma.log" 2>&1 &
  Rscript 06_fit_pirls_phi.R "$ex" moment  $MC_CORES > "run_${ex}_pirlsmoment.log"  2>&1 &
  Rscript 05_fit_lme4old.R   "$ex"         $MC_CORES > "run_${ex}_lme4old.log"      2>&1 &
  wait
  echo "--- $ex six-method extras done ---"
  tail -n 5 "run_${ex}_pirlsdigamma.log" "run_${ex}_pirlsmoment.log" "run_${ex}_lme4old.log"

  echo "### $ex : combined analysis"
  Rscript 07_analysis.R "$ex" > "run_${ex}_analysis6.log" 2>&1
  tail -n 25 "run_${ex}_analysis6.log"
  echo "=== $ex COMPLETE ==="
done

echo "############################################"
echo "### six-panel plot"
echo "############################################"
Rscript 08_summary_plots.R > run_08_summary_plots.log 2>&1
cat run_08_summary_plots.log

echo "############################################"
echo "### ALL SIX-METHOD EXTRAS COMPLETE"
echo "############################################"
