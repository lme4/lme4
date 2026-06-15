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

singularity exec \
    --bind "${RESULTS_DIR}:/results" \
    --env "REVDEP_LME4=${REVDEP_LME4}" \
    --env "_R_CHECK_FORCE_SUGGESTS_=false" \
    "${CONTAINER}" \
    Rscript /opt/revdep/check_one.R "${SLURM_ARRAY_TASK_ID}"
