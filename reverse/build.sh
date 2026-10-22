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
## --full        : also install Suggests of reverse dependencies recursively
##                 (larger image, ~6-7 GB vs ~2-3 GB, but more faithful to --as-cran)
## --no-bioc     : skip Bioconductor repositories entirely (faster, avoids flaky
##                 Bioc index fetches; image tagged with -no_bioc suffix)
## --docker-only : build the Docker image but skip the singularity build step
##                 (useful when Singularity/Apptainer is not installed locally;
##                 convert later with: singularity build lme4_revdep.sif docker-daemon://IMAGE)
## --refresh BASE : update the existing image BASE with Dockerfile.refresh
##                 instead of building from scratch (image tagged with
##                 -refreshYYYYMMDD suffix); use the same --full/--no-bioc
##                 flags as BASE was built with
WITH_SUGGESTS=false
WITH_BIOC=true
DOCKER_ONLY=false
REFRESH_BASE=""
REFRESH_TAG=""
SUGGESTS_TAG=""
BIOC_TAG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --full)        WITH_SUGGESTS=true; SUGGESTS_TAG="-full";  shift ;;
        --no-bioc)     WITH_BIOC=false;    BIOC_TAG="-no_bioc";  shift ;;
        --docker-only) DOCKER_ONLY=true;                          shift ;;
        --refresh)     REFRESH_BASE="$2"; REFRESH_TAG="-refresh$(date +%Y%m%d)"; shift 2 ;;
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
IMAGE="lme4-revdep:${VER1}_vs_${VER2}${SUGGESTS_TAG}${BIOC_TAG}${REFRESH_TAG}"
SIF="${SCRIPT_DIR}/lme4_revdep.sif"

echo "lme4 tarball 1 : $TGZ1  ($VER1)"
echo "lme4 tarball 2 : $TGZ2  ($VER2)"
echo "Docker image   : $IMAGE"
echo "Singularity    : $SIF"
echo "with Suggests  : $WITH_SUGGESTS"
echo "with Bioc      : $WITH_BIOC"
echo "docker only    : $DOCKER_ONLY"
echo "refresh base   : ${REFRESH_BASE:-(none, full build)}"

## Assemble a minimal build context in a temp directory so that stray
## lme4_*.tar.gz files that accumulate in reverse/ from local workflow runs
## are not picked up by the Dockerfile's COPY lme4_*.tar.gz line.
BUILDCTX="$(mktemp -d)"
trap 'rm -rf "$BUILDCTX"' EXIT

cp "$TGZ1" "$TGZ2" "$BUILDCTX/"
cp "$SCRIPT_DIR/Dockerfile" \
   "$SCRIPT_DIR/Dockerfile.refresh" \
   "$SCRIPT_DIR/setup_revdeps.R" \
   "$SCRIPT_DIR/check_one.R" \
   "$SCRIPT_DIR/slurm_submit.sh" \
   "$SCRIPT_DIR/slurm_job.sh" \
   "$BUILDCTX/"

## Build Docker image (setup_revdeps.R inside determines old vs new by version)
if [ -n "${REFRESH_BASE}" ]; then
    docker build \
        --build-arg BASE="${REFRESH_BASE}" \
        --build-arg REFRESH_STAMP="$(date -u +%Y%m%dT%H%M%S)" \
        --build-arg WITH_SUGGESTS="${WITH_SUGGESTS}" \
        --build-arg WITH_BIOC="${WITH_BIOC}" \
        -t "$IMAGE" -f "$BUILDCTX/Dockerfile.refresh" "$BUILDCTX"
else
    docker build \
        --build-arg WITH_SUGGESTS="${WITH_SUGGESTS}" \
        --build-arg WITH_BIOC="${WITH_BIOC}" \
        -t "$IMAGE" "$BUILDCTX"
fi

## Convert to Singularity .sif (skipped with --docker-only)
if [ "${DOCKER_ONLY}" = "true" ]; then
    echo ""
    echo "Docker image built: ${IMAGE}"
    echo "Skipping Singularity conversion (--docker-only).  To convert later:"
    echo "  singularity build ${SIF} docker-daemon://${IMAGE}"
    echo "  # or transfer with docker save and convert on the login node:"
    echo "  docker save -o lme4_revdep.tar ${IMAGE}"
else
    singularity build "$SIF" "docker-daemon://${IMAGE}"
    echo ""
    echo "Done.  Singularity image: $SIF"
    echo ""
    echo "To submit job arrays on Compute Canada:"
    echo "  bash ${SCRIPT_DIR}/slurm_submit.sh $SIF results_old old [--account=...]"
    echo "  bash ${SCRIPT_DIR}/slurm_submit.sh $SIF results_new new [--account=...]"
fi
