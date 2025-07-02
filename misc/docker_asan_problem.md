
* `r-devel-san`
* installation of dependent packages via `r2u`

## References

* https://dirk.eddelbuettel.com/code/sanitizers.html

## Run DE sanitizer example

```
docker pull rocker/r-devel-san
docker run --rm -ti rocker/r-devel-san bash
wget https://cran.r-project.org/src/contrib/sanitizers_0.1.1.tar.gz
RD CMD INSTALL sanitizers_0.1.1.tar.gz
RD -e 'library(sanitizers); print(stackAddressSanitize(42))'
```

Fails, properly.



```
git clone --depth=1 https://github.com/eddelbuettel/r2u.git
cp -r r2u/docker/jammy/build r-devel-san-r2u
cd r-devel-san-r2u 
sed -i -e 's#rocker/r-ubuntu:22.04#rocker/r-devel-san#' Dockerfile
## fart around with Dockerfile.  *Still* can't get r2u enabled on top of r-devel-san;
##   wait for lme4 dependent packages to get installed instead
docker build -t r-devel-san-r2u-lme4 .
```

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

=======
docker run --rm -ti r-devel-san-r2u-lme4 bash
wget https://cran.r-project.org/src/contrib/sanitizers_0.1.1.tar.gz
RD CMD INSTALL sanitizers_0.1.1.tar.gz
Rscript -e 'sanitizers::stackAddressSanitize()'
wget https://cran.r-project.org/src/contrib/lme4_1.1-35.3.tar.gz
RD CMD INSTALL lme4_1.1-35.3.tar.gz
```

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

https://superuser.com/questions/1130898/no-internet-connection-inside-docker-containers

```
docker network create --driver bridge common
docker run -it --network common ubuntu:latest bash
docker run --rm -ti --network common r-devel-san-r2u-lme4 bash
```

```r
install.packages("Rcpp")
remotes::install_github("RcppCore/Rcpp", ref = "bugfix/r-vector-start-specializations")
```

##

```
docker run --rm -ti --network common rocker/r-devel-san bash
apt update
apt upgrade -y
apt install -y cmake
echo "install.packages('nloptr')" | RD --slave
echo "install.packages('remotes')" | RD --slave
echo "remotes::install_github('astamm/nloptr')" | RD --slave
```


```r
library(nloptr)
 fr <- function(x) {   ## Rosenbrock Banana function
         100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
     }
(S <- bobyqa(c(0, 0, 0), fr, lower = c(0, 0, 0), upper = c(0.5, 0.5, 0.5)))
```

@edd suggests building `r-devel-san` on top of `rocker/r2u` (see `r-devel-san/Dockerfile`)

This is hard to rebuild ATM because my network is being crappy. Can I use `rocker/r-devel` and rebuild just the `lme4` package with the desired sanitizer flags enabled?

**locally**: add `-fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer` to `CXXFLAGS` in `.R/Makevars`

ugh, might need g++ 13 ("undefined symbol: __ubsan_vptr_type_cache")

```
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install gcc-13
```

Doesn't seem to be an easy way to fix this in Ubuntu 22.04 (and installing gcc-13 risks messing up the

Try inside `r-devel-san-r2u-lme4` (or `rocker/r2u`) instead?


Having trouble overriding the built-in `-fsanitize` option (wanted to update this so that I could recompile just the package, not have to rebuild the whole system, which is problematic at the moment ...)

https://stackoverflow.com/questions/49217539/package-build-ignores-makevars-flags
https://cran.r-project.org/doc/FAQ/R-exts.html
https://cran.r-project.org/doc/manuals/R-admin.html#Customizing-package-compilation


```
docker run --rm -ti r-devel-san-r2u-lme4 bash
## https://stackoverflow.com/questions/49217539/package-build-ignores-makevars-flags
export newflags="CXXFLAGS=-fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer -std=gnu++17 -I"/usr/local/lib/R/include" -DNDEBUG  -I'/usr/local/lib/R/site-library/Rcpp/include' -I'/usr/local/lib/R/site-library/RcppEigen/include' -I'/usr/lib/R/library/Matrix/include' -I/usr/local/include   -DNDEBUG -DEIGEN_DONT_VECTORIZE -fpic  -g -O2 -Wall -pedantic -mtune=native"
mkdir -p .R
echo "$newflags" > .R/Makevars
RD CMD build --no-build-vignettes lme4
## MAKEFLAGS="$newflags" RD CMD INSTALL lme4_1.1-35.9000.tar.gz
R_MAKEVARS_USER=.R/Makevars RD CMD INSTALL lme4_1.1-35.9000.tar.gz
## still doesn't override
```

https://cran.r-project.org/doc/FAQ/R-admin.html#Customizing-package-compilation

Doesn't seem to override.

https://www.reddit.com/r/C_Programming/comments/yzuxzr/help_with_fsanitizeaddress_output/

```
export MV=`RD RHOME`/etc/Makevars.site
echo "CXXFLAGS=-fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer -std=gnu++17 -lasan" > $MV
```

maybe? (includes *both* `-fsanitize` statements ...)

Try `wch1/r-debug` build? Looks like it has the `address` as well as the `undefined` sanitizers ...

```
/usr/local/lib/R/site-library/00LOCK-lme4/00new/lme4/libs/lme4.so: undefined symbol: __asan_option_detect_stack_use_after_return
```

&(^*&*^(&*^&(*^(&^(^*(^&(*&^(&^Y!@$##%@#!!!

```
docker pull wch1/r-debug:latest
docker run --rm -ti wch1/r-debug:latest bash
apt update
apt upgrade -y
apt install -y cmake
## r86695
## ONLY sanitizing address, not undefined?
RDsan -e "install.packages('remotes')"
RDsan -e "install.packages('lme4')"  ## this is taking forever ... why?
## https://github.com/actions/runner-images/issues/9524
sysctl vm.mmap_rnd_bits=28

## with RDsan
install.packages("remotes")
remotes::install_github("astamm/nloptr")
```

but now we have https://github.com/wch/r-debug/issues/37 (can't install packages!)

How about using R-devel (which presumably doesn't have any `-fsanitize`) and using that in `PKG_CXXFLAGS`) ?

## trying again with r-hub checkers

non-interactive, but ... ?

```
cat ~/admin/github_token_rhub
```


set credentials wih `gitcreds::gitcreds_set()`, run `rhub::rhub_check()` with the `clang-asan` platform. 


## 
```
cd r-san-devel
docker build -t r-devel-san-r2u:latest .
docker run --rm -ti --network common r-devel-san-r2u:latest bash
RD -e "install.packages('lme4')"
wget https://cran.r-project.org/src/contrib/lme4_1.1-35.3.tar.gz
RD CMD INSTALL lme4_1.1-35.3.tar.gz
```

https://reside-ic.github.io/blog/debugging-and-fixing-crans-additional-checks-errors/

