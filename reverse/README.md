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
| `build.sh` | Convenience wrapper: copies both lme4 tarballs into build context, `docker build`, `singularity build` → single `lme4_revdep.sif` |
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
export NEW=2.0-2
export OLD=2.0-1
```

### 1. Build the Singularity image (locally, needs Docker + Singularity)

`build.sh` can be run from any directory — it uses its own location to find
the `reverse/` build context and always writes `lme4_revdep.sif` there.
Tarball arguments (if given) are resolved relative to the current working
directory; with no arguments it auto-detects the two most recent
`lme4_*.tar.gz` files in `lme4/` (the parent of `reverse/`).

```bash
# Build the dev lme4 tarball if needed
cd ../..
R CMD build lme4 --compact-vignettes=both          # produces lme4_NEW.tar.gz
cd lme4/reverse

# Download the previous CRAN release for comparison
wget https://cran.r-project.org/src/contrib/lme4_${OLD}.tar.gz

# Build one image containing both versions (run from lme4/reverse/)
bash build.sh lme4_${OLD}.tar.gz lme4_${NEW}.tar.gz   # -> lme4_revdep.sif
```

If Docker is available locally but Singularity is not, export and convert
on the login node instead:

```bash
docker save -o lme4_revdep.tar lme4-revdep:OLD_vs_NEW
# transfer to Compute Canada, then:
singularity build lme4_revdep.sif docker-archive:lme4_revdep.tar
```

### 2. Transfer the image to Compute Canada

**Option A: direct file transfer via scp**

```bash
scp lme4_revdep.sif username@cedar.computecanada.ca:~/revdep/
```

**Option B: via Docker Hub**

```bash
# Locally: tag and push to Docker Hub (replace 'myuser' with your Docker Hub username)
docker tag lme4-revdep:OLD_vs_NEW myuser/lme4-revdep:OLD_vs_NEW
docker push myuser/lme4-revdep:OLD_vs_NEW

# On the Compute Canada login node: pull and convert
singularity pull lme4_revdep.sif docker://myuser/lme4-revdep:OLD_vs_NEW
```

Note: Docker Hub images are public by default. If the image should be kept
private, either use a private Docker Hub repository or use Option A instead.

### 3. Submit checking job arrays

```bash
# On the Compute Canada login node:
cd ~/revdep
bash slurm_submit.sh lme4_revdep.sif results_old old --account=def-yourpi
bash slurm_submit.sh lme4_revdep.sif results_new new --account=def-yourpi
```

Each array runs one task per reverse dependency (up to 50 concurrently).
SLURM logs go to `results_*/slurm_JOBID_TASKID.{out,err}`.
Check results go to `results_*/rdepends_PKGNAME.Rcheck/`.

### 4. Compare results

Once both arrays complete, run `checkChanges.R` as usual:

```bash
R --vanilla -f checkChanges.R --args \
    --old=results_old --new=results_new
```
