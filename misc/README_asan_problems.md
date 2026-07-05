# `Downdated VtV is not positive definite`: investigation notes

(by Claude, 5 July 2026)

## Background

Upgrading to Rcpp 1.1.2 started producing intermittent `Downdated VtV is not
positive definite` errors from `glmer()` fits. The failures were flaky:
`R CMD BATCH --vanilla` on a test file would often succeed while `R CMD check`
(or the same file run under more memory pressure) would fail, and independently
an ASan run on a Mac M1 CRAN-check machine flagged a heap-use-after-free in
the same code path. This document summarizes what's been found and fixed, and
what is still open.

## Fix already applied (commit `48508a84`)

**"fixed heap-use-after-free bug in merPredD initialization"**

In `R/AllClass.R`, `merPredD$initializePtr()` built the C++ object with:

```r
Ptr <<- .Call(merPredDCreate, as(X, "matrix"), Lambdat, ...)
```

`X` is already coerced to class `"matrix"` in `initialize()` (`X <<- as(X,
"matrix")`) and stored as a locked field. But `initializePtr()` re-wrapped it
in `as(X, "matrix")` *inline*, as a fresh call argument. Under R's
reference-counting semantics this can produce a **new, unprotected temporary**
R object (distinct from the field), which is what actually got passed to the
C++ constructor. `merPredD`'s C++ side stores this as a raw
`Eigen::Map<MatrixXd>` (`d_X`), aliasing the SEXP's memory directly with no
`PROTECT` and no persistent SEXP reference of its own. Once the `.Call()`
returned, nothing on the R side referenced that temporary any more, so it
could be (and was) garbage collected while `d_X` still pointed at it —
a classic heap-use-after-free, matching the ASan report exactly (`X`
allocated via `Rf_allocMatrix`/`model.matrix()`, freed by `RunGenCollect`,
read in `merPredD::updateXwts()` / `predModule.cpp:221`).

Every sibling argument to `merPredDCreate` (`Lambdat`, `V`, `VtV`, `Zt`, ...)
was passed as a bare field reference; `X` was the sole outlier. Every other
`*_Create` call in `AllClass.R` (`lm_Create`, `glm_Create`, `nls_Create`,
`glmFamily_Create`, `golden_Create`, `NelderMead_Create`) follows the
bare-field pattern correctly.

**Fix:** pass the field directly instead of re-coercing it:

```r
Ptr <<- .Call(merPredDCreate, X, Lambdat, ...)
```

This fix is real and necessary — it resolved the originally-reported failure
(the `cbpp` `gm2` overdispersion example in `R CMD check`'s "checking
examples" step went from ERROR to OK) — but **it turned out not to be
sufficient**. A second, still-unresolved heap corruption bug remains (see
below).

A NEWS entry was added to `inst/NEWS.Rd` under `CHANGES IN VERSION 2.0-3` /
`BUG FIXES` describing this fix.

## The bug that remains

### Symptom

`Downdated VtV is not positive definite` (thrown at `predModule.cpp:296-297`,
`if (d_RX.info() != Eigen::Success) Rcpp::stop(...)`), occurring
intermittently in `glmer()` fits that have nothing numerically unusual about
them (ordinary, well-conditioned models — not just the deliberately singular
case in `test-isSingular.R`, see caveat below).

### Confirmed via AddressSanitizer

Building `lme4.so` with `-fsanitize=address -fno-omit-frame-pointer -g` (and
running with `LD_PRELOAD=$(gcc -print-file-name=libasan.so)
ASAN_OPTIONS=detect_leaks=0`) reliably catches a **heap-use-after-free**:

```
READ of size 8 ...
  in Eigen::internal::product_evaluator<
       Eigen::Product<Eigen::DiagonalWrapper<Eigen::Map<VectorXd>>,
                       Eigen::Map<MatrixXd>>>::coeff(...)
     .../Eigen/src/Core/ProductEvaluators.h:925
  in lme4::merPredD::updateXwts(...) predModule.cpp:221
     (d_V = d_Xwts.asDiagonal() * d_X;)
  in internal_glmerWrkIter  external.cpp:274
  in pwrssUpdate            external.cpp:330
  in glmerLaplace           external.cpp:386

freed by thread T0 here:
  ReleaseLargeFreeVectors  memory.c:1167
  RunGenCollect            memory.c:1951
  R_gc_internal            memory.c:3237
  (an ordinary, unrelated GC cycle — not gctorture)

previously allocated by thread T0 here:
  Rf_allocVector3  memory.c:2894
  Rf_allocMatrix   array.c:236
  modelmatrix      .../stats/src/model.c:676   <-- model.matrix()
```

So it's `X` again (the fixed-effect design matrix from `model.matrix()`),
freed by a completely ordinary GC cycle, read from `d_X` inside
`updateXwts()` — **after the fix in `48508a84` was already applied**.

### Ruled out: R-devel regression

Reproduced the identical failure on stable **R 4.6.1**, not just the R-devel
snapshot originally used (2026-06-29 r90199). Built a fully separate package
library (`~/R/lib-4.6.1`: Matrix, Rcpp, RcppEigen, MASS, Rdpack, boot, minqa,
nlme, nloptr, reformulas, lattice, testthat from CRAN, lme4 rebuilt from
source with `--preclean`) to avoid any R-devel-compiled object files leaking
in. 5 runs of `test-covariance_structures.R`: 1 clean failure with the
identical signature. So this is **not** a recent R-devel regression — it
reproduces on the same stable release most users run.

### What's been ruled out / established about the mechanism

- **Not tied to one specific model.** Reproduced via the `cbpp` `gm2`
  overdispersion example, the plain `Contraception` random-slopes model
  (`gm <- glmer(use ~ age + urban + (1 + urban | district), ...)`, no
  structured-covariance syntax involved), a deliberately-singular Poisson GLMM
  in `test-isSingular.R`, and a model in `test-predict.R`. It is a general,
  timing-dependent memory bug, not a property of any one dataset/formula.
- **Not about `us()`/`cs()`/`diag()` structured covariance.** Verified `gm`
  (plain syntax) alone, fit twice in the same session, reproduces the crash
  on the second fit — structured-covariance syntax is not required.
- **`X`'s R-level field is never reassigned.** Grepped and confirmed only one
  `X <<-` assignment exists (in `initialize()`). Added temporary address-
  tracking instrumentation (via `.Internal(inspect(X))`) to `initialize()`,
  `initializePtr()`, `ptr()`, and the R-level `updateXwts()` method: the
  address is stable and consistent across every call, for the entire
  lifetime of a fit.
- **The C++ `d_X.data()` pointer is stable too.** Temporary `Rprintf` probes
  in the `merPredD` constructor and in `updateXwts()` (C++) confirm `d_X`
  is bound once at construction and never changes — hundreds of subsequent
  calls (both via the R-level `updateXwts()` method and via the internal
  C++-only PIRLS loop used when `compDev=TRUE`, the default) all report the
  identical address, matching construction. This was checked with
  `glmerControl(compDev=FALSE)` too (forces the R-mediated
  `RglmerWrkIter`/`pp$updateXwts()` path instead of the tight C++
  `glmerLaplace`/`pwrssUpdate` loop) — same result, same bug.
- **Despite all of the above being internally consistent, the backing
  memory still gets freed** while the `merPredD` instance is demonstrably
  still in active use (sometimes within the very *first* `glmer()` call of a
  session — a second call isn't even required). I.e., `X` becomes
  GC-unreachable even though `pp$X` (the reference-class field holding the
  external pointer's companion data) appears, by every check above, to still
  be a live, correctly-bound reference.
- **One confirmed before/after data point** (from an `rr` reverse-execution
  watchpoint, see below): the corrupted memory held the double value `1.0`
  (`0x3FF0000000000000`, plausibly the model matrix's intercept column)
  immediately before being overwritten with an unrelated value
  (`5.527444190707525`) by some other, unrelated allocation reusing the
  freed page.

**The open question:** *why* does R's GC ever decide `X` is unreachable, when
by ordinary R semantics `pp$X` (a binding inside the `merPredD` reference-class
instance/environment, itself kept alive via `rho$pp` where `rho` is the
`devfun` closure's environment, kept alive by `devfun` being a live variable
on the optimizer's call stack) should protect it for as long as the fit is
running? This has not yet been answered.

### Two things that are *not* this bug (for future readers)

1. **`test-isSingular.R:99`** deliberately simulates data from
   `theta = c(0, 0, 0, 1.0, 0.5, 0.3)` (several variance/correlation
   components exactly zero — genuinely singular by construction) specifically
   to check that `isSingular()` flags it. Hitting the same error message here
   may just be the **zero-tolerance PD check** (`predModule.cpp:296-297`, no
   tolerance/regularization) legitimately tripping on a model that's right at
   a singular boundary — a real robustness gap, but conceptually distinct
   from the memory-corruption bug. It has not been isolated from the memory
   bug with certainty, since both produce an identical error message; it's
   flagged here so it isn't confused with the main investigation.
2. **Hardening the PD check is not a fix for the memory bug.** Widening the
   tolerance or falling back to a regularized decomposition would only
   suppress the visible symptom. The underlying heap-use-after-free would
   remain and would instead silently corrupt fixed-effect/deviance results in
   some fraction of fits with no error raised at all — arguably worse than
   the current loud failure. The two problems should be fixed independently.

## Where the bug currently reproduces

It requires substantial memory-allocation pressure — small, isolated repro
scripts essentially never trigger it, even across many repeated fits with
explicit `gc()` calls. Known reproduction points, in decreasing order of
reliability:

- **`R CMD check` on the whole package** (examples + tests + vignettes run in
  sequence) — was the original, most reliable trigger.
- **`tests/testthat/test_check("lme4")`** (the full test suite via
  `AAAtest-all.R`) — reproduces reliably, hitting different tests/models on
  different runs (seen failing at `test-covariance_structures.R:396`/`404`
  and separately at `test-predict.R:678`).
- **`tests/testthat/test-covariance_structures.R`** run alone via
  `testthat::test_file()` — reproduces roughly 1 run in 3-5, at the plain
  `gm <- glmer(use ~ age + urban + (1 + urban | district), data =
  Contraception, family = binomial)` model (line ~396-404) or its
  `us()`/`cs()`/`diag()` siblings.
- **`vignettes/glmer.Rnw`, "toenailplot" chunk** (~line 807-822): fitting
  `glmer(outcome ~ time*treatment + (1 | patientID), data = toenail,
  family = binomial(link = "logit"))`. This one is notable: it fails
  **reliably (3/3)** when the vignette is rendered via
  `tools::buildVignette()`/knitr, but **never** (0/3+) when the identical
  code — including every preceding chunk in the same script, tangled via
  `knitr::purl()` — is run as a plain `Rscript`. Something about knitr's
  rendering machinery (graphics device allocation overhead is the leading
  suspect, untested) supplies enough extra allocation churn to trigger it
  where a plain script doesn't.
- Under **ASan** (`-fsanitize=address`, `LD_PRELOAD=libasan.so`), the failure
  rate goes up noticeably (ASan's allocator perturbs timing/layout enough to
  make the race easier to hit), which is how the crash details above were
  captured.
- **`gctorture()`/`gctorture2()` were tried and found impractical**: even
  `gctorture2(step=1)` (GC on every allocation) on a single, minimal `glmer()`
  fit did not complete within ~10 minutes, and lighter settings
  (`gctorture2(step=50)`, or a single `gc(full=TRUE)` at the suspected
  vulnerable moment) did not reproduce the bug even when run against
  deliberately-unfixed code. This is consistent with the bug being a
  *delayed-reuse* UAF (the read only misbehaves once the freed block is
  actually overwritten by something else) rather than the "collected
  immediately after last reference" pattern `gctorture` is designed to catch.
  **Conclusion: don't try to build a `gctorture`-based regression test for
  this bug** — it's neither fast nor reliable enough.

## Diagnostic tooling notes (for whoever picks this up next)

Getting a genuine backward-in-time trace of *why* `X` became unreachable
required `rr` record-and-replay plus reverse `gdb` debugging. This
environment (an AMD Zen CPU inside what turned out to be a KVM guest,
despite not announcing itself as one) needed several fixes along the way,
documented here in case they're needed again:

1. **`rr` doesn't work at all on Zen CPUs out of the box** — recording fails
   immediately unless the hardware "SpecLockMap" optimization is disabled
   (`sudo apt-get install msr-tools`, `sudo modprobe msr`, then a `wrmsr`
   call per the instructions at
   <https://github.com/rr-debugger/rr/wiki/Zen>).
2. **Ubuntu 22.04's packaged `rr` (5.5.0) can't record at all in this VM**:
   recording fails with `[record_syscall] File [vvar_vclock] is outside
   known root` / `EIO`. This is a known, already-fixed upstream issue
   (rr PR #3910, merged Jan 2025, not yet backported to the Ubuntu package).
   Fix: install a newer `rr` release directly
   (`rr-5.9.0-Linux-x86_64.deb` from
   <https://github.com/rr-debugger/rr/releases>).
3. **`gdb` cannot `continue` past an `execve()` replay event** — this is
   a known, still-open `rr` limitation
   (<https://github.com/rr-debugger/rr/issues/2381>, "Support executing past
   execve with gdb"). The maintainer's documented workaround: find an event
   number past the last exec (`rr replay -g <event> -p <pid> -s 0`,
   increasing `<event>` until the target executable shown in the "Launch
   debugger with ..." banner is the real one, e.g. `.../bin/exec/R` rather
   than `/usr/bin/bash`), then connect a debugger to that already-positioned
   session instead of starting from the very beginning of the recording.
4. **Ubuntu's packaged `gdb` (12.1) segfaults internally** (inside
   `bpstat_stop_status`, `breakpoint.c:5990`) immediately upon resuming a
   `-g`-restarted session with a Python-defined breakpoint active. Building
   `gdb` 17.2 from source (`./configure --prefix=~/gdb-17.2
   --with-python=python3 && make -j4 && make install`; needed only `flex`
   beyond what was already installed) did **not** fix this by itself.
5. **The actual cause of the gdb segfault was our own script**: the
   `gdb.Breakpoint.stop()` callback was calling
   `gdb.execute("delete breakpoints")`, which deletes the *currently
   executing* breakpoint object from inside its own callback — a
   use-after-free inside GDB's own breakpoint-processing code
   (`b = bs->breakpoint_at; ...; ++(b->hit_count);` right after our callback
   returns, on a since-freed `b`). Fix: use `self.enabled = False` instead of
   deleting. Combined with gdb 17.2, this reaches the crash site reliably
   (100% success across repeated attempts).
6. **`reverse-continue` runs asynchronously** over `target extended-remote`
   — issuing it and then immediately calling other commands (`bt`, setting
   another watchpoint) fails with "Cannot execute this command while the
   (selected thread/target) is running", because reverse execution hasn't
   actually finished yet. `set mi-async off` (tried both from within the
   Python script and before connecting, in the driver `.gdb` file) did not
   change this. The correct fix is to structure the walk as event-driven:
   register a `gdb.events.stop` handler and only issue the next
   `watch`/`reverse-continue` once the stop event for the previous one has
   actually arrived; the *first* `reverse-continue` also needs to be
   deferred out of the initiating breakpoint's `stop()` callback via
   `gdb.post_event()`, or it's rejected outright ("Cannot execute this
   command while the selected thread is running").
7. **`rr`'s reverse execution is not true hardware reverse-stepping** — it
   locates the nearest earlier checkpoint and replays forward repeatedly to
   find the exact target point, which can take a genuinely long time
   (multiple minutes) for a recording as large as a full `R CMD check`/test
   suite run. Give it a long timeout; what looks like a hang may just be
   in-progress work.
8. **The whole setup is flaky at the tooling level, independent of any of
   the above.** With every fix above applied, the identical
   record/attach/breakpoint/reverse-continue sequence has gone from
   working (reaching the crash site and completing one real
   `reverse-continue`, with genuine before/after watchpoint data — see next
   section) to failing instantly (~1 second, not a timeout, no error
   traceback) on many consecutive attempts, with no code change in between.
   One contributing factor found and fixed: a *stale, forgotten background
   retry-loop process* left running from an earlier attempt was still alive
   and periodically calling `pkill -9 -f "rr replay"`, which will race with
   and kill any other in-progress attempt's server. **Always check `ps aux
   | grep "rr replay\|rr_retry_loop"` for leftover processes from previous
   attempts before troubleshooting further** — but that was not confirmed
   to be the sole explanation, since instant failures were also observed
   immediately after clearing all such stale processes.

The scripts (`rr_watch.py`, `rr_retry_loop.sh`,
`rr_retry_loop_reverse_walk.sh`) and a fuller usage writeup live in
`misc/asan_rr_debug/` in this repo (not the session's ephemeral `/tmp`
scratch directory — those files will not survive past this session). Start
there. As of this writing, the event-driven `reverse-continue` walk (item 6
above) reaches the crash site reliably but the full backward walk to the
root cause has not yet completed/succeeded end-to-end, for the flakiness
reasons in item 8.

## Repository state / cleanup needed

- `src/Makevars` currently has debug flags added for this investigation
  (`-g -O0` appended to `PKG_CXXFLAGS`, no ASan flags currently active) —
  **should be reverted to the original** (`PKG_CXXFLAGS = -DNDEBUG
  -DEIGEN_DONT_VECTORIZE`, no `-g -O0`) before any release/commit.
- The fix in `R/AllClass.R` (bare `X` in `initializePtr()`) and the
  corresponding `inst/NEWS.Rd` entry are already committed (`48508a84`) and
  should be kept.
- No fix for the remaining bug has been written yet; nothing else in the
  source tree should need reverting, but double-check `git diff` against
  `48508a84` before committing anything further, since various temporary
  `Rprintf`/`cat()` debug instrumentation was added and removed multiple
  times during this investigation.
