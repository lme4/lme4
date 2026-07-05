#!/bin/bash
# Same idea as rr_retry_loop.sh, but the success criterion is that the
# event-driven reverse-continue walk in rr_watch.py actually received at
# least one stop event (i.e. reverse-continue genuinely completed at least
# once), not just that the crash site was reached.
#
# See rr_retry_loop.sh for the environment variables this reads
# (TRACE, TARGET_PID, GOTO_EVENT, GDB, WORKDIR, ATTEMPTS); this script adds:
#   GDB_TIMEOUT seconds to allow each gdb session (default: 240). rr's
#               reverse-continue is not true hardware reverse execution --
#               it locates the nearest earlier checkpoint and replays
#               forward repeatedly to find the exact point, which can take
#               a genuinely long time on a big recording. Increase this if
#               you suspect it just needs more time.

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$HERE/work}"
GDB="${GDB:-gdb}"
ATTEMPTS="${ATTEMPTS:-8}"
GDB_TIMEOUT="${GDB_TIMEOUT:-240}"

: "${TRACE:?set TRACE to an rr trace name/path, e.g. Rscript-10}"
: "${TARGET_PID:?set TARGET_PID to the recorded pid to attach to}"
: "${GOTO_EVENT:?set GOTO_EVENT to an event number past all execve()s for TARGET_PID}"

mkdir -p "$WORKDIR"

for i in $(seq 1 "$ATTEMPTS"); do
  echo "=== attempt $i ==="
  LOGF="$WORKDIR/rr_server_rw_$i.log"
  INVF="$WORKDIR/rr_gdb_rw_$i.log"
  pkill -9 -f "rr replay" 2>/dev/null
  sleep 1
  setsid nohup rr replay -g "$GOTO_EVENT" -p "$TARGET_PID" -s 0 -k "$TRACE" \
    > "$LOGF" 2>&1 < /dev/null &
  until grep -q "Launch debugger" "$LOGF" 2>/dev/null; do sleep 1; done
  PORT=$(grep -oP "127.0.0.1:\K[0-9]+" "$LOGF" | tail -1)
  TARGET_EXE=$(grep -oP "target extended-remote 127.0.0.1:[0-9]+' '\K[^']+" "$LOGF" | tail -1)
  cat > "$WORKDIR/rr_driver_rw_$i.gdb" << GDBEOF
set sysroot /
set mi-async off
target extended-remote 127.0.0.1:$PORT
source $HERE/rr_watch.py
continue
GDBEOF
  echo "START $(date +%s)" >> "$INVF"
  timeout "$GDB_TIMEOUT" "$GDB" -l 10000 -x "$WORKDIR/rr_driver_rw_$i.gdb" "$TARGET_EXE" \
    >> "$INVF" 2>&1
  echo "END $(date +%s)" >> "$INVF"
  if grep -q "stop event received\|REVERSE WALK FINISHED\|Old value" "$INVF"; then
    echo ">>> GOT REAL PROGRESS on attempt $i -- see $INVF <<<"
    break
  else
    echo "attempt $i: no reverse-continue completion -- see $INVF"
  fi
done
echo "RETRY_LOOP_DONE"
