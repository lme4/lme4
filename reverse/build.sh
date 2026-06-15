#!/usr/bin/env bash
##
## build.sh -- build a single Docker image containing both lme4 versions
##             and convert it to a Singularity .sif file
##
## Usage:
##   bash build.sh lme4_OLD.tar.gz lme4_NEW.tar.gz
##
## If no arguments are given, the two most recent lme4_*.tar.gz files in
## the parent directory are used (sorted by modification time, oldest first).
##
## The old and new tarballs are sorted by version number inside the image
## (see setup_revdeps.R) so argument order doesn't matter.
##
## Produces:
##   lme4_revdep.sif  in the reverse/ directory
##
## Requires Docker on the local machine.  The .sif conversion also requires
## Singularity/Apptainer locally.  If only Docker is available, use
## 'docker save' to export and convert on a machine that has Singularity
## (see the Dockerfile header for details).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

## ---- Parse flags ------------------------------------------------------------
## --full  : also install Suggests of reverse dependencies recursively
##           (larger image, ~6-7 GB vs ~2-3 GB, but more faithful to --as-cran)
WITH_SUGGESTS=false
SUGGESTS_TAG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --full) WITH_SUGGESTS=true; SUGGESTS_TAG="-full"; shift ;;
        --) shift; break ;;
        -*) echo "Unknown flag: $1"; exit 1 ;;
        *)  break ;;
    esac
done

if [ $# -ge 2 ]; then
    TGZ1="$(realpath "$1")"
    TGZ2="$(realpath "$2")"
else
    mapfile -t FOUND < <(ls -t "${PARENT_DIR}"/lme4_*.tar.gz 2>/dev/null | head -2)
    [ "${#FOUND[@]}" -eq 2 ] \
        || { echo "Need two lme4_*.tar.gz files; pass them as arguments or place them in ${PARENT_DIR}"; exit 1; }
    TGZ1="${FOUND[1]}"   # older (ls -t puts newest first)
    TGZ2="${FOUND[0]}"
fi

VER1="$(basename "$TGZ1" .tar.gz | sed 's/^lme4_//')"
VER2="$(basename "$TGZ2" .tar.gz | sed 's/^lme4_//')"
IMAGE="lme4-revdep:${VER1}_vs_${VER2}${SUGGESTS_TAG}"
SIF="${SCRIPT_DIR}/lme4_revdep.sif"

echo "lme4 tarball 1 : $TGZ1  ($VER1)"
echo "lme4 tarball 2 : $TGZ2  ($VER2)"
echo "Docker image   : $IMAGE"
echo "Singularity    : $SIF"
echo "with Suggests  : $WITH_SUGGESTS"

## Assemble a minimal build context in a temp directory so that stray
## lme4_*.tar.gz files that accumulate in reverse/ from local workflow runs
## are not picked up by the Dockerfile's COPY lme4_*.tar.gz line.
BUILDCTX="$(mktemp -d)"
trap 'rm -rf "$BUILDCTX"' EXIT

cp "$TGZ1" "$TGZ2" "$BUILDCTX/"
cp "$SCRIPT_DIR/Dockerfile" \
   "$SCRIPT_DIR/setup_revdeps.R" \
   "$SCRIPT_DIR/check_one.R" \
   "$SCRIPT_DIR/slurm_submit.sh" \
   "$SCRIPT_DIR/slurm_job.sh" \
   "$BUILDCTX/"

## Build Docker image (setup_revdeps.R inside determines old vs new by version)
docker build --build-arg WITH_SUGGESTS="${WITH_SUGGESTS}" -t "$IMAGE" "$BUILDCTX"

## Convert to Singularity .sif
singularity build "$SIF" "docker-daemon://${IMAGE}"

echo ""
echo "Done.  Singularity image: $SIF"
echo ""
echo "To submit job arrays on Compute Canada:"
echo "  bash ${SCRIPT_DIR}/slurm_submit.sh $SIF results_old old [--account=...]"
echo "  bash ${SCRIPT_DIR}/slurm_submit.sh $SIF results_new new [--account=...]"
