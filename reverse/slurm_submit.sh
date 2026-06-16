#!/usr/bin/env bash
##
## slurm_submit.sh -- submit a SLURM job array for lme4 revdep checking
##
## Usage:
##   bash slurm_submit.sh CONTAINER RESULTS_DIR old|new [extra sbatch options]
##
## The third argument selects which lme4 version to test against:
##   old  -> /opt/revdep/Library_old/  (the previous CRAN release)
##   new  -> /opt/revdep/Library_new/  (the dev version under test)
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
## The %50 throttle on the array limits concurrent tasks to 50; adjust to
## suit your allocation's fair-share policy.

set -euo pipefail

MYACCOUNT="${MYACCOUNT:-def-bolker}"

CONTAINER="${1:?Usage: $0 CONTAINER RESULTS_DIR old|new [extra sbatch options]}"
RESULTS_DIR="${2:?Usage: $0 CONTAINER RESULTS_DIR old|new [extra sbatch options]}"
LME4_VER="${3:?Usage: $0 CONTAINER RESULTS_DIR old|new [extra sbatch options]}"
shift 3   # remaining args forwarded to sbatch (e.g. --account=, --partition=)

[[ "$LME4_VER" == "old" || "$LME4_VER" == "new" ]] \
    || { echo "third argument must be 'old' or 'new'"; exit 1; }

CONTAINER="$(realpath "$CONTAINER")"
mkdir -p "$RESULTS_DIR"
RESULTS_DIR="$(realpath "$RESULTS_DIR")"

## Count packages directly from the container image
## (pipe through host wc -l; redirection inside singularity exec needs sh -c)
N=$(singularity exec "$CONTAINER" cat /opt/revdep/pkgs_to_check.txt | wc -l)
N="${N//[[:space:]]/}"
echo "Container      : $CONTAINER"
echo "Results dir    : $RESULTS_DIR"
echo "lme4 version   : $LME4_VER"
echo "Packages       : $N"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

sbatch \
    --array="1-${N}%50" \
    --time=2:00:00 \
    --mem=4G \
    --cpus-per-task=1 \
    --account="${MYACCOUNT}" \
    --job-name="lme4_revdep_${LME4_VER}" \
    --output="${RESULTS_DIR}/slurm_%A_%a.out" \
    --error="${RESULTS_DIR}/slurm_%A_%a.err" \
    --export="ALL,CONTAINER=${CONTAINER},RESULTS_DIR=${RESULTS_DIR},REVDEP_LME4=${LME4_VER}" \
    "$@" \
    "${SCRIPT_DIR}/slurm_job.sh"

echo "Job array submitted (1-${N}, lme4=${LME4_VER})."
echo "Results will appear in ${RESULTS_DIR}/"
