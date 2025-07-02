---
title: "trying to make sense of LAPACK/troubleshoot CRAN problems"
---

See also https://stat.ethz.ch/pipermail/r-package-devel/2022q2/007936.html

Rocker notes [here](https://www.rocker-project.org/), [here](https://github.com/rocker-org/rocker/issues/134), [here](https://github.com/r-lib/remotes/issues/371)

It's useful to use the `rocker/verse` image (has more stuff we need) rather than the `r-base` image.

```
docker pull rocker/verse ## fails intermittently?
docker run -ti rocker/verse bash
```

Maybe we need the `r-devel` container instead, which is built on Debian testing/unstable rather than stable, so that we can get a sufficiently recent OpenBLAS version?

```
docker pull rocker/r-devel
docker run -ti rocker/r-devel bash
```

From bash shell within rocker:

```
apt-get update ## maybe unnecessary (but doesn't take long)
install2.r --deps TRUE lme4
```

Needed to set up network by hand to get connectivity for installing packages etc. ???

from [here](https://docs.docker.com/engine/reference/commandline/network_connect/) etc.:
```
docker network create my-net
docker run -it --network=my-net rocker/r-devel bash
```

test:

```
cat >test.R <<EOF
set.seed(101)
d <- data.frame(z=rnorm(200),
                f=factor(sample(1:10,200, replace=TRUE)))

library(lme4)
library(testthat)
fm1 <- lmer(z~ as.numeric(f) + 1|f, d)
## need to \-protect $ signs
fm1@optinfo\$derivs\$Hessian[2,2] <- NA
expect_warning(lme4:::checkConv(fm1@optinfo\$derivs,
                                coefs=c(1,1),
                                ctrl=lmerControl()\$checkConv,lbound=0),
               "Problem with Hessian check")
EOF


cat >check_lapack.R <<EOF
s <- sessionInfo()
cat(s\$BLAS, s\$LAPACK, sep = "\n")
EOF
```

```
R --vanilla <check_lapack.R
R CMD BATCH --vanilla test.R; cat test.Rout
```

Now investigate/try to recreate the problem. From Kurt Hornik ...

> The failing check runs are using an external LAPACK 3.10.0. To be even more precise, Brian does not get this on Fedora with a system reference (netlib) 3.10 LAPACK, whereas I do on Debian with OpenBLAS.  Hth.


Default versions (from `sessionInfo()`):

```
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
```

```
update-alternatives --query liblapack.so.3-x86_64-linux-gnu
```

Choices: `openblas-pthread/liblapack.so.3` (default), `lapack/liblapack.so.3`


Manual:

```
update-alternatives --config libblas.so.3-x86_64-linux-gnu
update-alternatives --config liblapack.so.3-x86_64-linux-gnu
```

Now `LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0`

Fails to throw errors.


## Next questions: 

- is this an LAPACK/BLAS interaction, i.e. do I have to mess around with the BLAS version as well? (Seems unlikely.)
- (How) can I get so.3.10.0 on this machine? Do I have to move to a different/later Ubuntu version? (This is version 20.04)


## farting around with library-switching on local machine

So far I am baffled. What configure flags do I need to build R-devel with in order to allow switchable BLAS/LAPACK via `update-alternatives` (or some other method)? Do I need `--with-lapack --with-blas` or should I *not* include those alternatives? At present `update-alternatives` gives me choices of Intel MKL, OpenBLAS, or vanilla `libblas.so.3`, but running `update-alternatives` does **not** appear to change the reported results of `sessionInfo()` ... ?????

Why does Rocker report separate entries for `BLAS`/`LAPACK` while my machine reports a single `BLAS/LAPACK` value ... ?? (Rocker R version is 4.1.3 vs. r-devel on my machine ...)

## to test full package

```
wget https://cran.r-project.org/src/contrib/lme4_1.1-29.tar.gz
R CMD check --as-cran lme4_1.1-29.tar.gz
```

## more minimal tests

```r
dd <- list(gradient = c(0.00132136676711525, 0.00898197413334856, 0
), Hessian = structure(c(195.253128051758, 962.483270645142,
0, 962.483270645142, NA, 0, 0, 0, 6542.44775390625), dim = c(3L,
3L)))
with(dd, solve(chol(Hessian),gradient))
```


## junk

### attempt to get all the bits we need for LaTeX

definitely a nuisance, haven't sorted this all out yet

```
install2.r tinytex
## apt install texlive-latex-extra ## fails 'unmet dependencies'
apt install texlive-font-utils  ## gets all of texlive-base ??
```

updating TeX repos (from [here](https://tug.org/texlive/upgrade.html))

```
wget https://mirror.ctan.org/systems/texlive/tlnet/update-tlmgr-latest.sh
chmod +x update-tlmgr-latest.sh
./update-tlmgr-latest.sh
tlmgr update --self --all
```

installing way too many LaTeX packages:

```
tlmgr install fancyvrb
tlmgr install xcolor
tlmgr install natbib
tlmgr install amsmath
tlmgr install framed
tlmgr install hyperref
tlmgr install iftex
tlmgr install epstopdf
tlmgr install pdftexcmds
tlmgr install infwarerr
```

### For r-base only

Not needed with `rocker/verse` (which uses binary package installs and has more stuff installed by default)

```
apt install cmake  ## for nloptr from source
apt install qpdf
apt install libgsl-dev ## for gsl from source -> semEff (lme4 dep)
```

