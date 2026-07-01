#!/usr/bin/env bash
# Diff the current (HEAD) version of an .Rnw vignette against an older commit.
#
# Usage: ./latexdiff_commit.sh <commit-hash> [filename.Rnw]
#   filename.Rnw defaults to glmer.Rnw
set -euo pipefail

if [ -z "${1:-}" ]; then
    echo "Usage: $0 <commit-hash> [filename.Rnw]" >&2
    exit 1
fi

COMMIT="$1"
FN="${2:-glmer.Rnw}"
BASE="$(basename "$FN" .Rnw)"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "Working tree has uncommitted changes; commit or stash them first." >&2
    exit 1
fi

cleanup() {
    git switch - >/dev/null 2>&1 || true
}
trap cleanup ERR

git checkout --detach "$COMMIT"
Rscript -e "knitr::knit(\"$FN\")"
cp "${BASE}.tex" "${BASE}_old.tex"

git switch -

Rscript -e "knitr::knit(\"$FN\")"

latexdiff "${BASE}_old.tex" "${BASE}.tex" --math-markup=whole > "${BASE}_diff.tex"

texi2dvi -p "${BASE}_diff.tex"
