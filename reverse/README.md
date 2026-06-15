# lme4 reverse dependency checking

This directory contains scripts for checking reverse dependencies of lme4,
both locally (original workflow) and on Compute Canada HPC via
Singularity/Apptainer and SLURM (Docker-based workflow).

## Files

| File | Role |
|---|---|
| `Dockerfile` | Starts from `rocker/r2u`, installs TeX/tools, runs `setup_revdeps.R` at build time |
| `setup_revdeps.R` | Build-time: discovers rev deps via CRAN+Bioc, downloads source tarballs, installs all transitive deps + dev lme4 |
| `check_one.R` | Per-SLURM-task: runs `R CMD check --as-cran` on one tarball; names output `rdepends_PKG.Rcheck` for compatibility with `checkChanges.R` |
| `slurm_submit.sh` | Queries the container for N packages, submits an `sbatch --array=1-N%50` job |
| `slurm_job.sh` | Individual SLURM task: loads Apptainer, bind-mounts the results dir, calls `check_one.R` |
| `build.sh` | Convenience wrapper: finds the tarball, `docker build`, `singularity build` → `.sif` |
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

### 1. Build the Singularity images (locally, needs Docker + Singularity)

```bash
# Build R CMD build output first if needed
cd ..
R CMD build lme4          # produces lme4_VERSION.tar.gz
cd reverse

# Build one image per lme4 version to compare
bash build.sh ../lme4_2.0-1.tar.gz   # -> lme4_revdep_2.0-1.sif
bash build.sh ../lme4_2.0-2.tar.gz   # -> lme4_revdep_2.0-2.sif
```

If Docker is available locally but Singularity is not, export and convert
on the login node instead:

```bash
docker save -o lme4_revdep_2.0-2.tar lme4-revdep:2.0-2
# transfer to Compute Canada, then:
singularity build lme4_revdep_2.0-2.sif docker-archive:lme4_revdep_2.0-2.tar
```

### 2. Transfer images to Compute Canada

```bash
scp lme4_revdep_2.0-1.sif lme4_revdep_2.0-2.sif \
    username@cedar.computecanada.ca:~/revdep/
```

### 3. Submit checking job arrays

```bash
# On the Compute Canada login node:
cd ~/revdep
bash slurm_submit.sh lme4_revdep_2.0-1.sif results_old --account=def-yourpi
bash slurm_submit.sh lme4_revdep_2.0-2.sif results_new --account=def-yourpi
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
