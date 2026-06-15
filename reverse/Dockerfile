## Dockerfile for lme4 reverse dependency checking on HPC via Singularity/SLURM.
##
## Build (run from the reverse/ directory, after placing lme4_*.tar.gz here):
##   docker build -t lme4-revdep:VERSION .
##
## Convert to Singularity/Apptainer on a machine that has Docker:
##   singularity build lme4_revdep_VERSION.sif docker-daemon://lme4-revdep:VERSION
##
## If Docker is not available on the HPC login node, save/transfer instead:
##   docker save -o lme4_revdep.tar lme4-revdep:VERSION
##   # transfer to Compute Canada, then:
##   singularity build lme4_revdep_VERSION.sif docker-archive:lme4_revdep.tar
##
## Then submit the checking job array with slurm_submit.sh.

FROM rocker/r2u:latest

ENV DEBIAN_FRONTEND=noninteractive

## Tools needed by R CMD check --as-cran:
##   qpdf / ghostscript : PDF size/validity checks
##   texlive-*          : LaTeX-based vignettes
##   pandoc             : Rmarkdown vignettes
## r2u's bspm layer handles additional system-library deps for R packages
## automatically when install.packages() is called in setup_revdeps.R.
RUN apt-get update && apt-get install -y --no-install-recommends \
        qpdf \
        ghostscript \
        pandoc \
        texlive-latex-base \
        texlive-latex-recommended \
        texlive-latex-extra \
        texlive-fonts-recommended \
        texlive-fonts-extra \
    && rm -rf /var/lib/apt/lists/*

## Working directory for all revdep machinery
RUN mkdir -p /opt/revdep/tarballs

## Scripts
COPY setup_revdeps.R check_one.R slurm_submit.sh slurm_job.sh /opt/revdep/

## The lme4 source tarball to test (build context must contain lme4_*.tar.gz)
COPY lme4_*.tar.gz /opt/revdep/

## Internet-connected phase (runs at build time, not at check time):
##   - discover all CRAN + Bioconductor reverse dependencies of lme4
##   - download their source tarballs into /opt/revdep/tarballs/
##   - install transitive hard dependencies (Depends/Imports/LinkingTo only,
##     not recursive Suggests -- that balloons to thousands of packages)
##   - install both lme4 versions into Library_old/ and Library_new/
## Clean the apt cache in the same layer to avoid keeping 6+ GB of .deb files.
RUN Rscript /opt/revdep/setup_revdeps.R && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

CMD ["bash"]
