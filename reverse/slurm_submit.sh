#!/usr/bin/env bash
##
## slurm_submit.sh -- submit a SLURM job array for lme4 revdep checking
##
## Usage:
##   bash slurm_submit.sh CONTAINER RESULTS_DIR old|new [extra sbatch options]
##   bash slurm_submit.sh CONTAINER RESULTS_DIR both [extra sbatch options]
##
## The third argument selects which lme4 version to test against:
##   old  -> /opt/revdep/Library_old/  (the previous CRAN release)
##   new  -> /opt/revdep/Library_new/  (the dev version under test)
##   both -> a single array (not two) where each task checks the package
##           against old and new sequentially. Halves the number of jobs
##           submitted at once -- use this if your account hits
##           AssocMaxSubmitJobLimit submitting old and new separately
##           (a --dependency=afterany: chain does NOT help here: SLURM
##           counts dependency-held array tasks against the submit limit
##           too). Results land in RESULTS_DIR/old/ and RESULTS_DIR/new/.
##
## Example (mirroring the temprun workflow with a single shared image):
##
##   bash build.sh lme4_2.0-1.tar.gz lme4_2.0-2.tar.gz  # -> lme4_revdep.sif
##
##   bash slurm_submit.sh lme4_revdep.sif results_old old --account=def-yourpi
##   bash slurm_submit.sh lme4_revdep.sif results_new new --account=def-yourpi
##
##   # After both arrays complete:
##   R --vanilla -f checkChanges.R --args --old=results_old --new=results_new
##
## Or, submitting half as many jobs (see 'both' above):
##
##   bash slurm_submit.sh lme4_revdep.sif results both --account=def-yourpi
##   R --vanilla -f checkChanges.R --args --old=results/old --new=results/new
##
## The %50 throttle on the array limits concurrent tasks to 50; adjust to
## suit your allocation's fair-share policy.

set -euo pipefail

MYACCOUNT="${MYACCOUNT:-def-bolker}"
MAILUSER="${MAILUSER:-bolker@mcmaster.ca}"

CONTAINER="${1:?Usage: $0 CONTAINER RESULTS_DIR old|new [extra sbatch options]}"
RESULTS_DIR="${2:?Usage: $0 CONTAINER RESULTS_DIR old|new [extra sbatch options]}"
LME4_VER="${3:?Usage: $0 CONTAINER RESULTS_DIR old|new [extra sbatch options]}"
shift 3   # remaining args forwarded to sbatch (e.g. --account=, --partition=)

[[ "$LME4_VER" == "old" || "$LME4_VER" == "new" || "$LME4_VER" == "both" ]] \
    || { echo "third argument must be 'old', 'new', or 'both'"; exit 1; }

CONTAINER="$(realpath "$CONTAINER")"
mkdir -p "$RESULTS_DIR"
RESULTS_DIR="$(realpath "$RESULTS_DIR")"

if [[ "$LME4_VER" == "both" ]]; then
    RESULTS_DIR_OLD="${RESULTS_DIR}/old"
    RESULTS_DIR_NEW="${RESULTS_DIR}/new"
    mkdir -p "$RESULTS_DIR_OLD" "$RESULTS_DIR_NEW"
fi

## Compute Canada uses Apptainer (apptainer/singularity module); adjust the
## module name below to match your cluster ("apptainer" or "singularity").
## Needed here (not just in slurm_job.sh) because we call singularity
## directly below, on the login node, to count packages.
command -v singularity >/dev/null 2>&1 || module load apptainer

## Count packages directly from the container image
## (pipe through host wc -l; redirection inside singularity exec needs sh -c)
N=$(singularity exec "$CONTAINER" cat /opt/revdep/pkgs_to_check.txt | wc -l)
N="${N//[[:space:]]/}"

## All informational output goes to stderr so that stdout carries only the
## job ID, allowing clean capture: JOBID=$(bash slurm_submit.sh ...)
echo "Container      : $CONTAINER"  >&2
echo "Results dir    : $RESULTS_DIR" >&2
echo "lme4 version   : $LME4_VER"   >&2
echo "Packages       : $N"           >&2

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

## 'both' mode runs two R CMD check invocations per task instead of one, so
## give it double the wall-clock budget.
if [[ "$LME4_VER" == "both" ]]; then
    JOBTIME=8:00:00
    EXPORT_VARS="ALL,CONTAINER=${CONTAINER},RESULTS_DIR_OLD=${RESULTS_DIR_OLD},RESULTS_DIR_NEW=${RESULTS_DIR_NEW},REVDEP_LME4=${LME4_VER},CHECK_ONE_R=${SCRIPT_DIR}/check_one.R"
else
    JOBTIME=4:00:00
    EXPORT_VARS="ALL,CONTAINER=${CONTAINER},RESULTS_DIR=${RESULTS_DIR},REVDEP_LME4=${LME4_VER},CHECK_ONE_R=${SCRIPT_DIR}/check_one.R"
fi

## --parsable outputs only the numeric job ID to stdout
JOBID=$(sbatch --parsable \
    --array="1-${N}%50" \
    --time="${JOBTIME}" \
    --mem=16G \
    --cpus-per-task=1 \
    --account="${MYACCOUNT}" \
    --mail-type=BEGIN,FAIL \
    --mail-user="${MAILUSER}" \
    --job-name="lme4_revdep_${LME4_VER}" \
    --output="${RESULTS_DIR}/slurm_%A_%a.out" \
    --error="${RESULTS_DIR}/slurm_%A_%a.err" \
    --export="${EXPORT_VARS}" \
    "$@" \
    "${SCRIPT_DIR}/slurm_job.sh")

echo "Job array submitted: ${JOBID} (1-${N}, lme4=${LME4_VER})." >&2
if [[ "$LME4_VER" == "both" ]]; then
    echo "Results will appear in ${RESULTS_DIR_OLD}/ and ${RESULTS_DIR_NEW}/" >&2
else
    echo "Results will appear in ${RESULTS_DIR}/"                      >&2
    echo "NOTE: --dependency=afterany does NOT avoid AssocMaxSubmitJobLimit"  >&2
    echo "(dependency-held array tasks still count against the submit quota)." >&2
    echo "Either wait for this array to fully drain before submitting the"     >&2
    echo "other version, or resubmit using 'both' mode instead of 'old'/'new'." >&2
fi
echo "${JOBID}"
