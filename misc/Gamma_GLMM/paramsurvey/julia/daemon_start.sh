#!/usr/bin/env bash
# Start a persistent Julia DaemonMode.jl server with this directory's
# project (MixedModels.jl pa/dispersion-again + deps) activated, so
# packages load/precompile once and stay warm across script runs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

PORT="${1:-3141}"
PIDFILE="daemon.pid"

if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "daemon already running (pid $(cat "$PIDFILE"))"
  exit 0
fi

nohup julia --project=. --startup-file=no -e "using DaemonMode; serve($PORT)" \
  > daemon.log 2>&1 &
echo $! > "$PIDFILE"
echo "started daemon on port $PORT (pid $!), log: daemon.log"
