#!/usr/bin/env bash
##
## slurm_submit.sh -- submit a SLURM job array for lme4 revdep checking
##
## Usage:
##   bash slurm_submit.sh CONTAINER RESULTS_DIR [extra sbatch options]
##
## Example (new version vs old version workflow, mirroring temprun):
##
##   # Build containers for old and new lme4 (see Dockerfile):
##   bash build.sh lme4_2.0-1.tar.gz   # -> lme4_revdep_2.0-1.sif
##   bash build.sh lme4_2.0-2.tar.gz   # -> lme4_revdep_2.0-2.sif
##
##   # Submit checking job arrays
##   bash slurm_submit.sh lme4_revdep_2.0-1.sif ~/revdep/results_old
##   bash slurm_submit.sh lme4_revdep_2.0-2.sif ~/revdep/results_new
##
##   # After both arrays complete, compare (see checkChanges.R):
##   R --vanilla -f checkChanges.R --args \
##       --old=~/revdep/results_old --new=~/revdep/results_new
##
## The %50 throttle on the array limits concurrent tasks to 50; adjust to
## suit your allocation's fair-share policy.

set -euo pipefail

CONTAINER="${1:?Usage: $0 CONTAINER RESULTS_DIR [extra sbatch options]}"
RESULTS_DIR="${2:?Usage: $0 CONTAINER RESULTS_DIR [extra sbatch options]}"
shift 2   # remaining args forwarded to sbatch (e.g. --account=, --partition=)

CONTAINER="$(realpath "$CONTAINER")"
mkdir -p "$RESULTS_DIR"
RESULTS_DIR="$(realpath "$RESULTS_DIR")"

## Count packages directly from the container image
N=$(singularity exec "$CONTAINER" wc -l < /opt/revdep/pkgs_to_check.txt)
N="${N//[[:space:]]/}"
echo "Container      : $CONTAINER"
echo "Results dir    : $RESULTS_DIR"
echo "Packages       : $N"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

sbatch \
    --array="1-${N}%50" \
    --time=2:00:00 \
    --mem=4G \
    --cpus-per-task=1 \
    --job-name=lme4_revdep \
    --output="${RESULTS_DIR}/slurm_%A_%a.out" \
    --error="${RESULTS_DIR}/slurm_%A_%a.err" \
    --export="ALL,CONTAINER=${CONTAINER},RESULTS_DIR=${RESULTS_DIR}" \
    "$@" \
    "${SCRIPT_DIR}/slurm_job.sh"

echo "Job array submitted (1-${N})."
echo "Results will appear in ${RESULTS_DIR}/"
