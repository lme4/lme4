#!/usr/bin/env bash
##
## grep_problems.sh -- scan a results_* directory of rdepends_PKG.Rcheck
## outputs for likely early-abort / missing-system-dependency failures
## (e.g. a forgotten 'module load gsl'), rather than legitimate R CMD
## check ERROR/WARNING/NOTE results.
##
## Usage:
##   bash grep_problems.sh RESULTS_DIR [extra grep -E pattern to OR in]
##
## Two passes:
##   1. grep 00install.out and 00check.log in each *.Rcheck dir for common
##      missing-system-library/header/configure-failure signatures.
##   2. flag any 00check.log under SHORT_LOG_THRESHOLD lines -- a real
##      R CMD check --as-cran run is on the order of 100+ lines even for a
##      trivial package (see reverse/v2.0.3_local_check/*.Rcheck/00check.log
##      for reference), so a much shorter log is a sign R CMD check bailed
##      out early rather than actually running the full check.

set -euo pipefail

RESULTS_DIR="${1:?Usage: $0 RESULTS_DIR [extra grep -E pattern]}"
EXTRA_PATTERN="${2:-}"

PATTERN='cannot find -l|\.h(pp)?: No such file or directory|configure: error|unable to load shared object|error while loading shared librar|ld: cannot find|ERROR: dependency|ERROR: compilation failed|ERROR: configuration failed'
[ -n "$EXTRA_PATTERN" ] && PATTERN="${PATTERN}|${EXTRA_PATTERN}"

SHORT_LOG_THRESHOLD=30   # lines; see note above

echo "=== Packages with missing-system-dependency-style errors ==="
found_any=0
for d in "$RESULTS_DIR"/rdepends_*.Rcheck; do
    [ -d "$d" ] || continue
    pkg="$(basename "$d" .Rcheck)"; pkg="${pkg#rdepends_}"
    for logfile in "$d/00install.out" "$d/00check.log"; do
        [ -f "$logfile" ] || continue
        hit=$(grep -E -i -m1 "$PATTERN" "$logfile" || true)
        if [ -n "$hit" ]; then
            found_any=1
            echo "[$pkg] $(basename "$logfile"): $hit"
        fi
    done
done
[ "$found_any" -eq 0 ] && echo "(none found)"

echo ""
echo "=== Packages with a suspiciously short 00check.log (< $SHORT_LOG_THRESHOLD lines) ==="
found_any=0
for d in "$RESULTS_DIR"/rdepends_*.Rcheck; do
    [ -d "$d" ] || continue
    pkg="$(basename "$d" .Rcheck)"; pkg="${pkg#rdepends_}"
    logfile="$d/00check.log"
    if [ ! -f "$logfile" ]; then
        found_any=1
        echo "[$pkg] MISSING 00check.log entirely"
        continue
    fi
    n=$(wc -l < "$logfile")
    if [ "$n" -lt "$SHORT_LOG_THRESHOLD" ]; then
        found_any=1
        echo "[$pkg] $n lines"
    fi
done
[ "$found_any" -eq 0 ] && echo "(none found)"

echo ""
echo "=== Unique 'configure: error:' messages across all logs (count, message) ==="
grep -h -oE '^configure: error:.*' \
    "$RESULTS_DIR"/rdepends_*.Rcheck/00install.out \
    "$RESULTS_DIR"/rdepends_*.Rcheck/00check.log 2>/dev/null \
    | sort | uniq -c | sort -rn \
    || echo "(none found)"
