#!/usr/bin/env bash
##
## build.sh -- build Docker image for a given lme4 tarball and convert to .sif
##
## Usage:
##   bash build.sh [lme4_VERSION.tar.gz]
##
## If no tarball is given, the most recent lme4_*.tar.gz in the parent
## directory (i.e. built by 'R CMD build') is used.
##
## Produces:
##   lme4_revdep_VERSION.sif  in the reverse/ directory
##
## Requires Docker on the local machine.  The .sif conversion also requires
## Singularity/Apptainer locally.  If only Docker is available, use
## 'docker save' to export and convert on a machine that has Singularity
## (see the Dockerfile header for details).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

## Locate the lme4 tarball
if [ $# -ge 1 ]; then
    TARBALL="$(realpath "$1")"
else
    TARBALL="$(ls -t "${PARENT_DIR}"/lme4_*.tar.gz 2>/dev/null | head -1)"
    [ -n "$TARBALL" ] || { echo "No lme4_*.tar.gz found in ${PARENT_DIR}"; exit 1; }
fi

PKG_VER="$(basename "$TARBALL" .tar.gz | sed 's/^lme4_//')"
IMAGE="lme4-revdep:${PKG_VER}"
SIF="${SCRIPT_DIR}/lme4_revdep_${PKG_VER}.sif"

echo "lme4 tarball : $TARBALL  (version $PKG_VER)"
echo "Docker image : $IMAGE"
echo "Singularity  : $SIF"

## Copy tarball into Docker build context (reverse/ dir)
cp -f "$TARBALL" "$SCRIPT_DIR/"

## Build Docker image
docker build -t "$IMAGE" "$SCRIPT_DIR"

## Convert to Singularity .sif
singularity build "$SIF" "docker-daemon://${IMAGE}"

echo ""
echo "Done.  Singularity image: $SIF"
echo ""
echo "To submit the checking job array on Compute Canada:"
echo "  bash ${SCRIPT_DIR}/slurm_submit.sh $SIF /path/to/results [--account=...]"
