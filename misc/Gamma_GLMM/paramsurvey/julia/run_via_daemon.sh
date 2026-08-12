#!/usr/bin/env bash
# Run a Julia script against the warm DaemonMode server started by
# daemon_start.sh, instead of paying full package-load/precompile cost
# for every invocation.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

PORT="${DAEMON_PORT:-3141}"
SCRIPT="$1"
shift || true

julia --project=. --startup-file=no -e "using DaemonMode; runargs($PORT)" "$SCRIPT" "$@"
