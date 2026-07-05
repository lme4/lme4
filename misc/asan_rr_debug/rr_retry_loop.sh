#!/bin/bash
# Repeatedly attach gdb to an rr replay session (past the execve boundary,
# see README.md) until CrashBreak in rr_watch.py reports hitting the
# "Downdated VtV" throw site, or ATTEMPTS is exhausted.
#
# Configure via environment variables (see README.md for how to determine
# TRACE / TARGET_PID / GOTO_EVENT for a given recording):
#   TRACE       rr trace name/path, e.g. Rscript-10 (relative to
#               ~/.local/share/rr/) or an absolute path.  Required.
#   TARGET_PID  pid (as originally recorded) of the process to attach to.
#               Required.
#   GOTO_EVENT  event number to jump to with `rr replay -g`, chosen to land
#               after all execve()s have completed for TARGET_PID. Required.
#   GDB         path to the gdb binary to use (default: gdb from PATH; on
#               the machine this was developed on, needed a from-source
#               build at ~/gdb-17.2/bin/gdb -- see README.md).
#   WORKDIR     directory for logs/generated driver files (default: ./work
#               next to this script).
#   ATTEMPTS    number of attempts before giving up (default: 10).
#   GDB_TIMEOUT seconds to allow each gdb session (default: 300).

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$HERE/work}"
GDB="${GDB:-gdb}"
ATTEMPTS="${ATTEMPTS:-10}"
GDB_TIMEOUT="${GDB_TIMEOUT:-300}"

: "${TRACE:?set TRACE to an rr trace name/path, e.g. Rscript-10}"
: "${TARGET_PID:?set TARGET_PID to the recorded pid to attach to}"
: "${GOTO_EVENT:?set GOTO_EVENT to an event number past all execve()s for TARGET_PID}"

mkdir -p "$WORKDIR"

for i in $(seq 1 "$ATTEMPTS"); do
  echo "=== attempt $i ==="
  LOGF="$WORKDIR/rr_server_$i.log"
  INVF="$WORKDIR/rr_gdb_$i.log"
  pkill -9 -f "rr replay" 2>/dev/null
  sleep 1
  setsid nohup rr replay -g "$GOTO_EVENT" -p "$TARGET_PID" -s 0 "$TRACE" \
    > "$LOGF" 2>&1 < /dev/null &
  until grep -q "Launch debugger" "$LOGF" 2>/dev/null; do sleep 1; done
  PORT=$(grep -oP "127.0.0.1:\K[0-9]+" "$LOGF" | tail -1)
  TARGET_EXE=$(grep -oP "target extended-remote 127.0.0.1:[0-9]+' '\K[^']+" "$LOGF" | tail -1)
  cat > "$WORKDIR/rr_driver_$i.gdb" << GDBEOF
set sysroot /
set mi-async off
target extended-remote 127.0.0.1:$PORT
source $HERE/rr_watch.py
continue
GDBEOF
  timeout "$GDB_TIMEOUT" "$GDB" -l 10000 -x "$WORKDIR/rr_driver_$i.gdb" "$TARGET_EXE" \
    > "$INVF" 2>&1
  if grep -q "CRASH SITE HIT" "$INVF"; then
    echo ">>> SUCCESS on attempt $i -- see $INVF <<<"
    break
  else
    echo "attempt $i failed (crash or no hit) -- see $INVF"
  fi
done
echo "RETRY_LOOP_DONE"
