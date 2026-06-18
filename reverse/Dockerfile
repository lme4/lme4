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
##   texinfo            : provides texi2dvi (required by some vignettes)
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
        texlive-lang-french \
        lmodern \
        texinfo \
    && rm -rf /var/lib/apt/lists/*

## Working directory for all revdep machinery
RUN mkdir -p /opt/revdep/tarballs

## Scripts
COPY setup_revdeps.R check_one.R slurm_submit.sh slurm_job.sh /opt/revdep/

## The lme4 source tarball to test (build context must contain lme4_*.tar.gz)
COPY lme4_*.tar.gz /opt/revdep/

## WITH_SUGGESTS=true  : include Suggests in recursive dep resolution (~6-7 GB)
## WITH_SUGGESTS=false : Depends/Imports/LinkingTo only          (~2-3 GB)
## WITH_BIOC=false     : skip Bioconductor repos (avoids flaky index fetches)
ARG WITH_SUGGESTS=false
ARG WITH_BIOC=true

## Internet-connected phase (runs at build time, not at check time):
##   - discover all CRAN + Bioconductor reverse dependencies of lme4
##   - download their source tarballs into /opt/revdep/tarballs/
##   - install transitive dependencies (depth controlled by WITH_SUGGESTS)
##   - install both lme4 versions into Library_old/ and Library_new/
## Clean the apt cache in the same layer to avoid keeping 6+ GB of .deb files.
RUN WITH_SUGGESTS=${WITH_SUGGESTS} WITH_BIOC=${WITH_BIOC} Rscript /opt/revdep/setup_revdeps.R && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

## Extra R packages needed to build certain vignettes but not declared as
## reverse-dependency dependencies (add to this list as needed).
RUN Rscript -e 'install.packages(c("fpc", "lmPerm"))'

CMD ["bash"]
