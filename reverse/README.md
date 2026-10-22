# lme4 reverse dependency checking

This directory contains scripts for checking reverse dependencies of lme4,
both locally (original workflow) and on Compute Canada HPC via
Singularity/Apptainer and SLURM (Docker-based workflow).

## Files

| File | Role |
|---|---|
| `Dockerfile` | Starts from `rocker/r2u`, installs TeX/tools, runs `setup_revdeps.R` at build time |
| `setup_revdeps.R` | Build-time: discovers rev deps via CRAN+Bioc, downloads source tarballs, installs all transitive deps, installs both lme4 versions into separate `Library_old/` and `Library_new/` dirs |
| `check_one.R` | Per-SLURM-task: prepends the appropriate versioned library, runs `R CMD check --as-cran` on one tarball; names output `rdepends_PKG.Rcheck` for compatibility with `checkChanges.R` |
| `slurm_submit.sh` | Queries the container for N packages, accepts `old`/`new` switch, submits an `sbatch --array=1-N%50` job |
| `slurm_job.sh` | Individual SLURM task: loads Apptainer, bind-mounts the results dir, passes `REVDEP_LME4` env var, calls `check_one.R` |
| `build.sh` | Convenience wrapper: copies both lme4 tarballs into build context, `docker build`, `singularity build` → single `lme4_revdep.sif`; `--refresh BASE` updates an existing image instead |
| `Dockerfile.refresh` | One layer on top of an existing image: runs `setup_revdeps.R` with `REFRESH=true` to bring revdeps and installed packages up to date |
| `Dockerfile.patch` | One layer on top of an existing image: reinstalls only the dev lme4 (plus a few fixed-up packages) |
| `checkReverse.R` | Original local workflow: installs deps and runs `R CMD check` via `tools::check_packages_in_dir` |
| `checkChanges.R` | Compares two sets of check results (old vs new lme4) using `tools::check_packages_in_dir_changes` |

## Local workflow

Mirrors the steps in `temprun`:

```bash
oldver=2.0-1
newver=2.0-2
R --vanilla -f checkReverse.R --args --jobs=4 lme4_${oldver}.tar.gz
R --vanilla -f checkReverse.R --args --jobs=4 lme4_${newver}.tar.gz
R --vanilla -f checkChanges.R --args \
    --old=lme4_${oldver}.tar.gz.reverse \
    --new=lme4_${newver}.tar.gz.reverse
```

## HPC workflow (Compute Canada / SLURM)

Compute Canada nodes may not have internet access, so the workflow is split
into a build phase (internet-connected, runs locally) and a check phase
(offline, runs on the cluster).

Both lme4 versions share a single container image: all reverse dependency
tarballs and their transitive dependencies are stored once, and each version
of lme4 is pre-installed into its own library directory (`Library_old/` and
`Library_new/`) inside the image.  Two job arrays then run against the same
image, selecting the appropriate library via the `REVDEP_LME4` environment
variable.

### 0. define versions

(example, useful if we want to cut and paste code below)
```bash
export NEW=2.1-0
export OLD=2.0-6
export CCUSER=bolker
export CCACCOUNT=def-bolker
export DOCKERUSER=bbolker
export TAG=${OLD}_vs_${NEW}-$(date +%Y%m%d)
export SIF=/scratch/$USER/lme4_revdep_${TAG}.sif
```

`TAG` is the name the image is pushed under (step 2). Builds of the same
`OLD`/`NEW` pair recur (the dev lme4 and the reverse dependencies keep
changing), so a dated tag keeps a new push from silently replacing the
Docker Hub image an earlier run used, and makes it obvious which image a
given `.sif` came from. Set `TAG` (and hence `SIF`) once per build/push and
reuse the same values on the cluster; `date` there will not necessarily
give the same day.

**Always refer to the `.sif` by this full absolute path (`$SIF`), never a bare
relative `lme4_revdep.sif`.** This `reverse/` directory accumulates images
from *previous, unrelated* version comparisons across releases (e.g. an old
image sitting there from checking a completely different `OLD`/`NEW` pair
months earlier) -- a bare relative filename silently resolves to whichever
`lme4_revdep.sif` happens to be in your current directory, and there is no
warning if that's stale. Also keep the `.sif` out of `$HOME` (small quota,
see the cache note below) -- `/scratch/$USER/` is a good default location.

### 1. Build the Singularity image (locally, needs Docker + Singularity)

`build.sh` can be run from any directory — it uses its own location to find
the `reverse/` build context and always writes `lme4_revdep.sif` there.
Tarball arguments (if given) are resolved relative to the current working
directory; with no arguments it auto-detects the two most recent
`lme4_*.tar.gz` files in `lme4/` (the parent of `reverse/`).

```bash
# Build the dev lme4 tarball if needed
cd ../..
R CMD build lme4 --compact-vignettes=both          # produces lme4_${NEW}.tar.gz
mv lme4_${NEW}.tar.gz lme4/reverse
cd lme4/reverse

# Download the previous CRAN release for comparison
wget https://cran.r-project.org/src/contrib/lme4_${OLD}.tar.gz

# Build one image containing both versions (run from lme4/reverse/)
bash build.sh lme4_${OLD}.tar.gz lme4_${NEW}.tar.gz   # -> lme4_revdep.sif
```

#### Build flags

| Flag | Effect | Image tag suffix | Approx. compressed size | Approx. build time |
|---|---|---|---|---|
| *(none)* | Deps/Imports/LinkingTo only; includes Bioconductor repos | *(none)* | ~6 GB | ~48 min |
| `--full` | Also resolves Suggests recursively | `-full` | ~25 GB | — |
| `--no-bioc` | Skip Bioconductor repositories (avoids flaky index fetches) | `-no_bioc` | ~4 GB | — |

Flags can be combined: `bash build.sh --full --no-bioc ...` produces a tag like
`lme4-revdep:${OLD}_vs_${NEW}-full-no_bioc`.

The default build includes Bioconductor but does not recurse into Suggests —
this is the recommended starting point. Use `--no-bioc` if the Bioconductor
index is consistently unreachable during build. Use `--full` only if you need
Suggests installed for the checks to be meaningful (significantly larger image
and longer build time).

If Docker is available locally but Singularity is not, export and convert
on the login node instead
(in this chunk and others below, container names should also include any tags, e.g. -full, -no_bioc):

```bash
docker save -o lme4_revdep.tar lme4-revdep:${OLD}_vs_${NEW}
# transfer to Compute Canada, then:
singularity build lme4_revdep.sif docker-archive:lme4_revdep.tar
```

#### Diagnosing build problems

`setup_revdeps.R` writes two files into the image at build time:

- `/opt/revdep/setup_summary.txt` — high-level summary (counts, lme4 versions, download failures)
- `/opt/revdep/install_failures.txt` — packages that failed to install (absent if none)

Inspect them without starting a full container:

```bash
docker run --rm lme4-revdep:${OLD}_vs_${NEW} cat /opt/revdep/setup_summary.txt
docker run --rm lme4-revdep:${OLD}_vs_${NEW} cat /opt/revdep/install_failures.txt
```

Or on the cluster from the `.sif` file:

```bash
singularity exec lme4_revdep.sif cat /opt/revdep/setup_summary.txt
singularity exec lme4_revdep.sif cat /opt/revdep/install_failures.txt
```

If packages failed to install, the compilation error usually names the missing
system library.  To reproduce it interactively:

```bash
docker run --rm -it lme4-revdep:${OLD}_vs_${NEW} bash
# inside the container:
R -e 'install.packages("PackageName")'
```

The error output (e.g. `hdf5.h: No such file or directory`) maps directly to
an apt package (e.g. `libhdf5-dev`).  Add it to the `apt-get install` block in
the Dockerfile and rebuild.

#### Updating an existing image (patch layer)

If the version pair is unchanged but the dev lme4 has moved on (or a system
library/package is missing), `Dockerfile.patch` adds one layer on top of the
existing image instead of a full rebuild: it reinstalls the new tarball into
`Library_new/` (leaving `Library_old/` alone), installs extra system
libraries/packages, and appends the lme4 git SHA to
`/opt/revdep/setup_summary.txt`. The revdep list and all other packages stay
as they were at the original build -- refresh the image (next section) or
do a full rebuild if those need updating.

Build the tarball from a clean checkout (e.g. a `git worktree`) so untracked
files in the working tree don't end up in it, then use a minimal build
context (as `build.sh` does) rather than `reverse/` itself:

```bash
export TAG=${OLD}_vs_${NEW}-patch$(date +%Y%m%d)   # overrides step 0's TAG
CTX=$(mktemp -d)
cp Dockerfile.patch lme4_${NEW}.tar.gz "$CTX"/
docker build --build-arg BASE=lme4-revdep:${OLD}_vs_${NEW} \
    --build-arg NEW_TGZ=lme4_${NEW}.tar.gz \
    --build-arg LME4_SHA=$(git rev-parse --short HEAD) \
    -t lme4-revdep:${TAG} -f "$CTX"/Dockerfile.patch "$CTX"
docker run --rm lme4-revdep:${TAG} cat /opt/revdep/setup_summary.txt
```

Push/pull it as in step 2.

#### Refreshing an existing image (updated revdeps and packages)

If the image is out of date with CRAN/Bioconductor (new reverse
dependencies, new versions of existing ones or of their dependencies) but
R itself hasn't changed minor version, `build.sh --refresh BASE` adds one
layer (`Dockerfile.refresh`) on top of an existing image instead of
rebuilding from scratch. It reruns `setup_revdeps.R` with `REFRESH=true`,
which

- re-discovers the reverse dependencies, keeps tarballs whose version is
  still current, downloads new/updated ones, and deletes the rest;
- holds `r-base-core` and runs `apt-get upgrade`, which updates only the
  r2u binary packages that changed;
- updates outdated packages in the source-built library
  (`/usr/local/lib/R/site-library`, mostly Bioconductor) from source, and
  reinstalls source-built packages whose `LinkingTo` dependencies (Rcpp,
  RcppEigen, TMB, ...) were upgraded, since they were compiled against the
  old headers;
- installs newly needed packages with the usual two passes;
- reinstalls both lme4 versions from the tarballs given to `build.sh`.

```bash
bash build.sh --docker-only --refresh lme4-revdep:${OLD}_vs_${NEW} \
    lme4_${OLD}.tar.gz lme4_${NEW}.tar.gz
# -> lme4-revdep:${OLD}_vs_${NEW}-refreshYYYYMMDD
```

Use the same `--full`/`--no-bioc` flags the base image was built with.
`setup_summary.txt` gets a `refreshed` line with counts of tarballs
downloaded/removed and packages upgraded/relinked.

For comparison (2026-10-22, refreshing an image built 2026-08-20): the
refresh took 19 min and added a 2.23 GB (uncompressed) layer, vs. ~50 min
for a full rebuild; it downloaded 141 tarballs, removed 110, upgraded 309 R
packages via apt and 3 from source, relinked 3, and installed 103 new ones.
Revdep tarballs and installed package versions then matched a full build
from the day before, apart from CRAN/r2u updates in between.

Caveats:

- Packages are never removed, so a refreshed image also contains packages
  that are no longer needed (36 in the example above). If one is in a
  reverse dependency's Suggests, that package's checks may run code a fresh
  build would skip; this affects old and new alike, so it shouldn't change
  the comparison.
- Every upgraded package is stored again in the new layer (the old copy
  stays in the base layers), so the image grows with each refresh. Do a
  full rebuild from time to time, and whenever R's minor version changes.
- The refresh reruns on every build (`REFRESH_STAMP` defeats the Docker
  cache), so it always picks up the current state of the repositories.
- Array task indices map to different packages than in the base image
  (the package list changes), so start the check run with a fresh results
  directory.

### 2. Transfer the image to Compute Canada

**Option A: direct file transfer via scp**

```bash
scp lme4_revdep.sif username@cedar.computecanada.ca:$SIF
```

**Option B: via Docker Hub**

```bash
# Locally: log in, tag, and push (the local image name is whatever
# build.sh printed, e.g. lme4-revdep:${OLD}_vs_${NEW}, plus any -full etc.)
docker login
docker tag lme4-revdep:${OLD}_vs_${NEW} ${DOCKERUSER}/lme4-revdep:${TAG}
docker push ${DOCKERUSER}/lme4-revdep:${TAG}

# On the Compute Canada login node (with OLD/NEW/TAG/SIF set as above):
# pull and convert, then check it's the image you meant to pull
module load apptainer/1.4.5
mkdir -p "$(dirname "$SIF")"
singularity pull "$SIF" docker://${DOCKERUSER}/lme4-revdep:${TAG}
singularity exec "$SIF" cat /opt/revdep/setup_summary.txt
```

`singularity pull` refuses to overwrite an existing file; `--force` makes it
do so, but prefer a new `$SIF` name per `TAG` (see step 0).

`singularity pull` caches downloaded Docker layers and conversion intermediates under
`$HOME/.apptainer` (or `$APPTAINER_CACHEDIR` if set), which can easily reach tens of GB for
an image this size -- `$HOME`'s quota on Compute Canada is small, and this cache is a common
way to blow it (symptom: obscure "Disk quota exceeded" errors on unrelated tiny file writes,
e.g. `.git/FETCH_HEAD`, once the filesystem is full or an inode limit is hit). Run
`diskusage_report` to check. Once the `.sif` is built you no longer need the cache, so it's
safe to just clear it:

```bash
apptainer cache clean   # or: rm -rf ~/.apptainer/cache
```

To keep the cache around (e.g. for future rebuilds) without it eating `$HOME`, relocate it and
point `APPTAINER_CACHEDIR` at the new location (add the `export` to `~/.bashrc` too, or it
reverts to `$HOME/.apptainer` on the next pull):

```bash
mkdir -p ~/project/$CCUSER/apptainer_cache
mv ~/.apptainer/cache/* ~/project/$CCUSER/apptainer_cache/
export APPTAINER_CACHEDIR=~/project/$CCUSER/apptainer_cache
```

### 3. Submit checking job arrays

A full run's `R CMD check` output (per package: install log, check log, example
`.Rout`s, unpacked sources, vignette builds, ...) generates a large number of
small files -- for ~600 reverse dependencies times two lme4 versions, easily
tens of thousands. `/project`'s per-file (inode) quota is comparatively small
and this can exhaust it outright even while byte usage still looks fine
(`diskusage_report` shows files and bytes as separate quotas -- check both).
Point the results-dir argument at `/scratch/$USER` instead of the `lme4/reverse`
checkout under `/project`, since `/scratch/$USER` has a much larger file-count
quota and the checked-out package tarballs don't need `/project`'s
persistence:

Recommended: submit a single combined array with `both` mode, where each
task checks its package against old and new sequentially. This halves the
number of jobs submitted at once compared to two separate arrays -- worth
doing by default, since many clusters' per-user `AssocMaxSubmitJobLimit` is
easy to hit otherwise. Note that a SLURM `--dependency` chain does **not**
help with that limit: dependency-held array tasks still count against the
submit quota the moment they're submitted, so staggering old/new via
`--dependency=afterany:` does not avoid it, it just delays when you hit it.

```bash
# On the Compute Canada login node:
cd ~/project/${CCUSER}/lme4/reverse
bash slurm_submit.sh "$SIF" /scratch/$USER/results both --account=${CCACCOUNT}
```

Results land in `/scratch/$USER/results/old/` and `/scratch/$USER/results/new/`.

If you'd rather run old and new as two independent arrays (e.g. to compare
timing, or because `both` mode's doubled per-task runtime doesn't fit your
`--time` budget), submit them separately:

```bash
bash slurm_submit.sh "$SIF" /scratch/$USER/results_old old --account=${CCACCOUNT}
bash slurm_submit.sh "$SIF" /scratch/$USER/results_new new --account=${CCACCOUNT}
```

If this hits the submit-job limit, the only real fix is to wait for the
first array to fully drain (`squeue -j JOBID -r` returns empty) before
submitting the second -- not `--dependency`, per above.

Each array runs one task per reverse dependency (up to 50 concurrently).
SLURM logs go to `results*/slurm_JOBID_TASKID.{out,err}`.
Check results go to `results*/rdepends_PKGNAME.Rcheck/` (nested under
`old/`/`new/` for `both` mode).

### 4. Compare results

Once the array (or both arrays) complete, run `checkChanges.R` as usual:

```bash
# both mode:
R --vanilla -f checkChanges.R --args \
    --old=/scratch/$USER/results/old --new=/scratch/$USER/results/new

# separate old/new arrays:
R --vanilla -f checkChanges.R --args \
    --old=/scratch/$USER/results_old --new=/scratch/$USER/results_new
```
