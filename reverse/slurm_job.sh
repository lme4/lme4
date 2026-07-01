#!/usr/bin/env bash
##
## slurm_job.sh -- single task in the lme4 revdep check job array
##
## Do not invoke directly; submitted by slurm_submit.sh.
## Expects these variables to be set (via --export in sbatch):
##   CONTAINER    : absolute path to the Singularity/Apptainer .sif file
##   RESULTS_DIR  : absolute path to the directory for .Rcheck outputs
##   REVDEP_LME4  : "old" or "new" -- selects which lme4 library dir to use
##   CHECK_ONE_R  : absolute path to check_one.R on the host filesystem
## SLURM_ARRAY_TASK_ID is set automatically by SLURM.
##
## Compute Canada uses Apptainer (apptainer/singularity module); adjust the
## module name below to match your cluster ("apptainer" or "singularity").

module load apptainer

## CHECK_ONE_R is passed via --export in slurm_submit.sh (absolute path on the
## host filesystem).  We bind-mount it over the baked-in copy so that script
## changes take effect without rebuilding the .sif image.
## Note: BASH_SOURCE[0] cannot be used here because SLURM copies the job
## script to a temporary spool directory, so check_one.R would not be found
## alongside it.

singularity exec \
    --bind "${RESULTS_DIR}:/results" \
    --bind "${CHECK_ONE_R}:/opt/revdep/check_one.R" \
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
