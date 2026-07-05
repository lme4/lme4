# rr + gdb reverse-execution scripts for the `Downdated VtV` UAF

Companion scripts to `../README_asan_problems.md`. These automate attaching
`gdb` to an `rr` replay session, past the point where the `merPredD::X`
use-after-free throws `Downdated VtV is not positive definite`, and (in the
`_reverse_walk` variant) walking backward through every write to the
corrupted memory via `reverse-continue` + a hardware watchpoint.

## Prerequisites (as of this writing, on the dev machine)

- **`rr` >= 5.9.0.** The Ubuntu-packaged 5.5.0 cannot even record in this
  environment (fails with `[vvar_vclock] is outside known root` / `EIO` —
  see main README's tooling notes). Grab a newer `.deb` directly from
  <https://github.com/rr-debugger/rr/releases>.
- **AMD Zen CPU workaround applied**, if applicable
  (<https://github.com/rr-debugger/rr/wiki/Zen>) — `rr record` fails
  immediately otherwise.
- **`gdb` >= ~17** (or otherwise a build where the segfault described below
  is fixed). Ubuntu 22.04's packaged `gdb` (12.1) segfaults internally
  (`bpstat_stop_status`, `breakpoint.c:5990`) as soon as a Python-defined
  breakpoint's `stop()` callback disables/deletes itself and a `-g`-restarted
  session resumes. Building from source
  (`./configure --prefix=~/gdb-17.2 --with-python=python3 && make -j4 &&
  make install`) worked; point `GDB=~/gdb-17.2/bin/gdb` at it if your system
  gdb is old.
- lme4 built with debug symbols (`-g` in `src/Makevars`'s `PKG_CXXFLAGS`;
  `-O0` recommended for cleaner line-level stepping, not required).
  **Do not build with ASan for this** — `LD_PRELOAD`ing `libasan.so` into
  gdb's own process hangs gdb; run the plain (non-ASan) debug build under
  `rr` instead. ASan is for the separate, non-`rr` reproduction described in
  the main README.

## Step 1: record a crashing run

The bug is intermittent, so record something likely to trigger it (the full
test suite has the best hit rate observed so far) and keep recording until
one fails:

```sh
cd /path/to/lme4/tests
Rscript -e 'library(testthat); library(lme4); test_check("lme4")' &
# or record test_file("test-covariance_structures.R") for a faster/lighter
# repro; roughly 1-in-3-to-5 runs fail.
```

under `rr record -n <that command>`. Once you get a trace whose output
contains `Downdated VtV is not positive definite`, note its name via
`rr ls` (traces live in `~/.local/share/rr/`; `rr ls | tail` or the
`latest-trace` symlink will point at the most recent one).

Confirm it replays deterministically:

```sh
rr replay -a <trace-name>   # should reproduce the same failure, no debugger
```

## Step 2: find TARGET_PID and GOTO_EVENT

`gdb` cannot `continue` past an `execve()` replay event — this is a known,
still-open `rr` limitation (rr issue #2381). Attach only *after* the
recorded process has already exec'd into the real target (e.g. R itself),
not into an intermediate shell:

```sh
rr ps <trace-name>          # lists pids/commands in the recording;
                             # TARGET_PID is normally the top-level one
                             # (the Rscript/R process), which keeps the same
                             # pid across its own internal exec() calls
```

Then binary-search for an event number past the last `execve`:

```sh
rr replay -g <EVENT> -p <TARGET_PID> -s 0 <trace-name>
```

Watch the "Launch debugger with ..." banner it prints: increase `<EVENT>`
until the executable named there is the real target (e.g.
`/usr/local/lib/R/bin/exec/R`) rather than an intermediate shell. `50000`
worked for a full-test-suite recording on the dev machine; this will vary by
trace. Kill that instance (`Ctrl-C` / `pkill -9 -f "rr replay"`) once you've
found a working event number — the scripts below manage the server
themselves.

## Step 3: run the scripts

```sh
export TRACE=<trace-name>            # or an absolute path
export TARGET_PID=<pid from rr ps>
export GOTO_EVENT=<event found above>
export GDB=~/gdb-17.2/bin/gdb        # if system gdb is too old

./rr_retry_loop.sh                   # just confirms CrashBreak fires
./rr_retry_loop_reverse_walk.sh      # attempts the full reverse walk
```

Logs and generated `.gdb` driver files land in `./work/` (gitignored).
Both scripts retry because **this has been observed to be flaky at the
tooling level** — the identical setup has, in one session, gone from
reaching the crash site and completing one real `reverse-continue` (with a
genuine before/after value at the watched address) to failing instantly on
many consecutive attempts with no code change (confirmed not explained by
a stale leftover background process holding a competing `rr replay` server
open — that was found and fixed once, but instant failures persisted
afterward too). If a run fails instantly (~1 second, not a timeout), that's
this flakiness, not a real result; just retry. If it fails by *hanging* for
the full timeout, that may be `reverse-continue` legitimately still
searching — see the checkpoint theory below.

### Diagnosed one layer deeper: the `bc` packet gets no response

Adding `set debug remote 1` to the driver (prints every gdb-remote-protocol
packet) shows exactly where an instant failure happens. Everything up to
and including arming the watchpoint succeeds normally
(`$Z2,<addr>,8#..` / `Packet received: OK`, "write-watchpoint is
supported"), then:

```
[remote] Sending packet: $bc#c5      <- 'bc' = reverse-continue
[remote] wait: enter
[remote] wait: exit                  <- no "Packet received" logged at all
Detaching from program: ...
❌️ Cannot execute this command while the target is running.
```

GDB's wait-for-stop-reply loop exits with **no packet received in between**
— rr's stub isn't sending back a proper stop-reply (or is dropping the
connection) in response to `bc`, and GDB reacts by detaching.

**Leading theory (untested):** `rr replay -g <event>` may fast-forward to
that event by a path that doesn't lay down the periodic checkpoints that
`reverse-continue` needs to search backward from — i.e. there may be
nothing behind the current position for `bc` to reverse *into*. If so, the
one earlier success was likely a fluke of some checkpoint existing to have
been created incidentally. This would mean **the whole `-g <event> -p <pid>`
workaround for the execve-follow limitation is fundamentally at odds with
doing a reverse walk** — it lets you attach past the exec boundary, but
maybe leaves you without the checkpoint history reverse execution depends
on.

**Tried and ruled out:** recorded a much smaller targeted trace
(`test_file("test-covariance_structures.R")` alone instead of the full test
suite — a `~30x` smaller `GOTO_EVENT`, 5000 vs. 50000, landing in the real
`R` binary just as quickly relative to the trace) and reran the reverse
walk against it. **Identical result: 8/8 instant failures**, same "no
packet received after `bc`" pattern. So it is *not* simply about how far
`-g` fast-forwards, or about needing "more room" for checkpoints to
accumulate before the target event — a 10x-smaller jump made no
difference at all. The checkpoint-scarcity theory above is probably wrong,
or at least incomplete.

**Not yet tried:**

1. See if `rr replay` has an option to explicitly force a checkpoint before
   attaching (rather than however `-g` positions things internally).
2. Try reverse-continue against a session that was attached *without* `-g`
   at all — i.e. reached the target point via genuinely stepping/continuing
   forward from the start of a recording (tolerating however `gdb` behaves
   on the handful of `execve`s along the way, e.g. by scripting past them
   with repeated top-level `continue`s rather than relying on `-g`/`-p` to
   skip them). This is the one variable that hasn't been isolated yet: every
   attempt so far, successful or not, used the `-g <event> -p <pid>`
   attach path. If the *one* prior success and all subsequent failures both
   used that path, the flakiness is happening somewhere else entirely (rr
   itself, this VM's environment, or system load at the time) — but this
   hasn't been tested against a non-`-g` baseline to know whether `-g`
   itself is even a contributing factor at all.
3. It may be worth searching/filing an rr GitHub issue with the exact
   protocol trace above (`bc` sent, no reply, gdb detaches) — this is
   specific and reproducible enough to be a good bug report, and the
   `rr` maintainers would know immediately whether `-g`-restarted sessions
   are expected to support `reverse-continue` at all.

If you get further with any of these, please update this file.

## What `rr_watch.py` does

- Sets a breakpoint on `predModule.cpp:297` (the `Rcpp::stop("Downdated VtV
  is not positive definite")` throw site inside `merPredD::updateDecomp`).
- On hit, reads `this->d_X.m_data` (the raw pointer backing the Eigen `Map`
  aliasing R's `X` matrix memory) via `gdb.parse_and_eval` — `.data()`
  itself often can't be called on optimized/inlined builds, hence the
  `m_data` fallback.
- **Disables itself (`self.enabled = False`) rather than deleting itself.**
  Deleting the currently-executing breakpoint from inside its own `stop()`
  callback was the actual root cause of the gdb-12.1 segfault mentioned
  above — GDB's C++ side still holds a reference to the breakpoint object
  after the Python callback returns and dereferences it, so deleting it
  first is a use-after-free inside GDB itself.
- Defers the actual reverse-continue walk via `gdb.post_event()` — issuing
  `reverse-continue` directly from within the breakpoint's own `stop()`
  callback is rejected ("Cannot execute this command while the selected
  thread is running"); it must run from a normal top-level context.
- Walks backward via a `gdb.events.stop`-driven state machine (not a plain
  sequential loop): `reverse-continue` runs asynchronously over
  `target extended-remote`, so the code must wait for the actual stop event
  rather than assuming `gdb.execute("reverse-continue")` blocks until done.
  (`set mi-async off` was tried as a fix for the same symptom and did not
  help.)
- Creates the watchpoint via the `watch` CLI command
  (`gdb.execute("watch ...")`), not the `gdb.Breakpoint(...,
  type=gdb.BP_WATCHPOINT)` Python constructor — the latter appears
  correlated with (though not conclusively proven to cause) the instant,
  no-error "Detaching from program" failures described above. Prefer the
  string-command form until/unless that's investigated further.

## What we already learned before the tooling got in the way

One earlier (non-reproducible-on-demand) successful run captured a real
watchpoint hit: the corrupted memory held the double value `1.0`
(`0x3FF0000000000000`, plausibly the model matrix's intercept column)
immediately before being overwritten with an unrelated value
(`5.527444190707525`) by some other allocation reusing the freed page. That
confirms the corruption mechanism but not yet *why* R's GC ever considered
`X` unreachable in the first place — that's the open question these scripts
are trying to answer.
