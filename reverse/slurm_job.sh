#!/usr/bin/env bash
##
## slurm_job.sh -- single task in the lme4 revdep check job array
##
## Do not invoke directly; submitted by slurm_submit.sh.
## Expects these variables to be set (via --export in sbatch):
##   CONTAINER    : absolute path to the Singularity/Apptainer .sif file
##   RESULTS_DIR  : absolute path to the directory for .Rcheck outputs
##   REVDEP_LME4  : "old" or "new" -- selects which lme4 library dir to use
## SLURM_ARRAY_TASK_ID is set automatically by SLURM.
##
## Compute Canada uses Apptainer (apptainer/singularity module); adjust the
## module name below to match your cluster ("apptainer" or "singularity").

module load apptainer

## Resolve the directory containing this script so we can bind-mount
## check_one.R from the host over the baked-in copy in the container.
## This means script changes take effect without rebuilding the .sif image.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

singularity exec \
    --bind "${RESULTS_DIR}:/results" \
    --bind "${SCRIPT_DIR}/check_one.R:/opt/revdep/check_one.R" \
    --env "REVDEP_LME4=${REVDEP_LME4}" \
    --env "_R_CHECK_FORCE_SUGGESTS_=false" \
    --env "_R_CHECK_CRAN_INCOMING_=false" \
    --env "R_PROFILE=/dev/null" \
    "${CONTAINER}" \
    Rscript /opt/revdep/check_one.R "${SLURM_ARRAY_TASK_ID}"
## R_PROFILE=/dev/null skips /usr/lib/R/etc/Rprofile.site, which calls
## bspm::enable().  bspm requires D-Bus (unavailable in Singularity on HPC)
## and emits a noisy warning; since all packages are pre-installed in the
## image and compute nodes have no internet, bspm serves no purpose here.
##
## _R_CHECK_CRAN_INCOMING_=false disables the entire CRAN incoming feasibility
## section of R CMD check --as-cran.  This covers both the remote feasibility
## sub-check (version already on CRAN, etc.) and the package-dependency check
## that fetches CRAN/Bioconductor PACKAGES indices -- both of which produce
## spurious network-access warnings on offline compute nodes.
