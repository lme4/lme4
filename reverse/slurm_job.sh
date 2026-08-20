#!/usr/bin/env bash
##
## slurm_job.sh -- single task in the lme4 revdep check job array
##
## Do not invoke directly; submitted by slurm_submit.sh.
## Expects these variables to be set (via --export in sbatch):
##   CONTAINER    : absolute path to the Singularity/Apptainer .sif file
##   REVDEP_LME4  : "old", "new", or "both" -- which lme4 library dir(s) to use
##   CHECK_ONE_R  : absolute path to check_one.R on the host filesystem
## For REVDEP_LME4 in {old, new}:
##   RESULTS_DIR  : absolute path to the directory for .Rcheck outputs
## For REVDEP_LME4 == both (halves the number of array tasks submitted,
## for clusters with a low AssocMaxSubmitJobLimit -- each task runs both
## versions sequentially instead of two separate array tasks doing one
## version each):
##   RESULTS_DIR_OLD, RESULTS_DIR_NEW : absolute paths for each version's
##                                      .Rcheck outputs
## SLURM_ARRAY_TASK_ID is set automatically by SLURM.
##
## Compute Canada uses Apptainer (apptainer/singularity module); adjust the
## module name below to match your cluster ("apptainer" or "singularity").
## gsl is a system dependency for many packages (packages are installed
## from source during the checking process).

module load apptainer
module load gsl

## CHECK_ONE_R is passed via --export in slurm_submit.sh (absolute path on the
## host filesystem).  We bind-mount it over the baked-in copy so that script
## changes take effect without rebuilding the .sif image.
## Note: BASH_SOURCE[0] cannot be used here because SLURM copies the job
## script to a temporary spool directory, so check_one.R would not be found
## alongside it.

run_one () {
    local ver="$1" results="$2"
    ## --no-home: without this, Singularity auto-binds the real $HOME (and
    ## sets HOME to it inside the container), so R's R_LIBS_USER picks up
    ## the host's personal library (e.g. ~/R/x86_64-pc-linux-gnu-library/*)
    ## ahead of the container's own r2u-installed packages in .libPaths().
    ## Those host .so files were built outside the container (different R
    ## build/glibc) and fail to load -- this shadows perfectly good
    ## in-container installs (stringi, mvtnorm, TMB, ragg, ...) with
    ## broken host ones. The explicit --bind mounts below are unaffected.
    singularity exec --no-home \
        --bind "${results}:/results" \
        --bind "${CHECK_ONE_R}:/opt/revdep/check_one.R" \
        --env "REVDEP_LME4=${ver}" \
        --env "_R_CHECK_FORCE_SUGGESTS_=false" \
        --env "_R_CHECK_CRAN_INCOMING_=false" \
        --env "R_PROFILE=/dev/null" \
        "${CONTAINER}" \
        Rscript /opt/revdep/check_one.R "${SLURM_ARRAY_TASK_ID}"
}

if [[ "${REVDEP_LME4}" == "both" ]]; then
    run_one old "${RESULTS_DIR_OLD}"
    run_one new "${RESULTS_DIR_NEW}"
else
    run_one "${REVDEP_LME4}" "${RESULTS_DIR}"
fi
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
