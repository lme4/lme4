```
docker pull ghcr.io/r-hub/containers/clang-asan
docker run --rm -ti ghcr.io/r-hub/containers/clang-asan
```

```
apt update
apt -y upgrade
apt -y install cmake apt-file wget libgsl-dev libpng-dev libjpeg-dev
apt-file update
Rscript -e "install.packages('remotes', repos = getOption('repos')[[2]])"

Rscript -e "remotes::install_github('cran4linux/rspm')"

## use bspm instead?? ugh.
## stolen/modified from Dockerfile for r/bspm, but rhub containers are different from rocker containers
## export RHOME=`Rscript -e "cat(Sys.getenv('R_HOME'))"`
## rm -f /etc/apt/sources.list.d/c2d4u*list \
##         && apt update -qq \ 
##         && Rscript -e "install.packages('bspm', repos = getOption('repos')[[2]])" \ 
##         && echo "suppressMessages(bspm::enable())" >> $RHOME/etc/Rprofile.site \ 
##         && echo "options(bspm.version.check=FALSE)" >> $RHOME/etc/Rprofile.site \ 
##         && echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/90local-no-recommends \ 
##         && chgrp -R 1000 $RHOME/library 

## Rscript -e "options(repos='https://packagemanager.rstudio.com/all/__linux__/focal/latest'); install.packages('lme4', dependencies = TRUE)"
Rscript -e "install.packages('lme4', dependencies = TRUE, repos = getOption('repos')[[2]])"

wget https://cran.r-project.org/src/contrib/lme4_1.1-35.3.tar.gz
R CMD check --as-cran lme4_1.1-35.3.tar.gz
```

## creating new docker file

```
cd R/pkgs/lme4/misc/r-devel-san-r2u
docker build -t r-devel-san-r2u-lme4 .
docker image ls
docker run --rm -ti r-devel-san-r2u-lme4 bash
RD CMD check --as-cran lme4_1.1-35.3.tar.gz 
```

```
https://dirk.eddelbuettel.com/code/sanitizers.html


```
## Install key and setup r2u repo, also set 'pin preference'
        && wget -q -O - https://r2u.stat.illinois.edu/ubuntu/dirk_eddelbuettel_pubkey.asc \
                | tee -a /etc/apt/trusted.gpg.d/dirk_eddelbuettel_pubkey.asc \
        && echo "deb [arch=amd64 signed-by=/etc/apt/trusted.gpg.d/dirk_eddelbuettel_pubkey.asc] https://r2u.stat.illinois.edu/ubuntu jammy main" \
                > /etc/apt/sources.list.d/r2u.list \
        && echo "Package: *" > /etc/apt/preferences.d/99r2u \
        && echo "Pin: release o=CRAN-Apt Project" >> /etc/apt/preferences.d/99r2u \
        && echo "Pin: release l=CRAN-Apt Packages" >> /etc/apt/preferences.d/99r2u \
        && echo "Pin-Priority: 700"  >> /etc/apt/preferences.d/99r2u \
## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
        && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8
```
