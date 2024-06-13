## Emacs, make this -*- mode: sh; -*-

## start with the Docker 'base R' Debian-based image + r2u
FROM rocker/r2u:latest

## This handle reaches Carl and Dirk
MAINTAINER "Carl Boettiger and Dirk Eddelbuettel" rocker-maintainers@eddelbuettel.com

## Remain current
RUN apt-get update -qq \
	&& apt-get dist-upgrade -y

## From the Build-Depends of the Debian R package, plus subversion
RUN apt-get update -qq \
    ## removed -t unstable in modifying from Debian to Ubuntu base
	&& apt-get install -y --no-install-recommends \
		bash-completion \
		bison \
		debhelper \
		default-jdk \
		g++ \
		gcc \
		gfortran \
		groff-base \
		libblas-dev \
		libbz2-dev \
		libcairo2-dev \
		libcurl4-openssl-dev \
		libjpeg-dev \
		liblapack-dev \
		liblzma-dev \
		libncurses-dev \
		libpango1.0-dev \
		libpcre2-dev \
		libpng-dev \
		libreadline-dev \
		libtiff5-dev \
		libx11-dev \
		libxt-dev \
		mpack \
		subversion \
		tcl8.6-dev \
		texinfo \
		texlive-base \
		texlive-extra-utils \
		texlive-fonts-extra \
		texlive-fonts-recommended \
		texlive-plain-generic \
		texlive-latex-base \
		texlive-latex-extra \
		texlive-latex-recommended \
		tk8.6-dev \
		valgrind \
		x11proto-core-dev \
		xauth \
		xdg-utils \
		xfonts-base \
		xvfb \
		zlib1g-dev 

## Check out R-devel
RUN cd /tmp \
	&& svn co https://svn.r-project.org/R/trunk R-devel 

## Build and install according the standard 'recipe' I emailed/posted years ago
## Updated compiler flags to match https://www.stats.ox.ac.uk/pub/bdr/memtests/README.txt
RUN cd /tmp/R-devel \
	&& R_PAPERSIZE=letter \
	   R_BATCHSAVE="--no-save --no-restore" \
	   R_BROWSER=xdg-open \
	   PAGER=/usr/bin/pager \
	   PERL=/usr/bin/perl \
	   R_UNZIPCMD=/usr/bin/unzip \
	   R_ZIPCMD=/usr/bin/zip \
	   R_PRINTCMD=/usr/bin/lpr \
	   LIBnn=lib \
	   AWK=/usr/bin/awk \
	   CC="gcc -std=gnu99 -fsanitize=undefined -fno-omit-frame-pointer" \
	   FFLAGS="-g -O2 -mtune=native" \
	   FCFLAGS="-g -O2 -mtune=native" \
	   CXX="g++-13 -fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer" \
	   CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native" \
	   CXX11FLAGS="-g -O2 -Wall -pedantic -mtune=native" \
	   CXX14FLAGS="-g -O2 -Wall -pedantic -mtune=native" \
	   CXX17FLAGS="-g -O2 -Wall -pedantic -mtune=native" \
	   MAIN_LDFLAGS="-fsanitize=undefined -pthread" \
	   F77="gfortran -fsanitize=undefined" \
	   FC="gfortran -fsanitize=undefined" \
	   FCFLAGS="-g -O2 -mtune=native" \
	   FFLAGS="-g -O2 -mtune=native" \
	   ./configure --enable-R-shlib \
	       --without-blas \
	       --without-lapack \
	       --with-readline \
	       --without-recommended-packages \
	       --program-suffix=dev \
	       --disable-openmp \
	&& make \
	&& make install \
	&& make clean

## Set Renviron to get libs from base R install
RUN echo "R_LIBS=\${R_LIBS-'/usr/local/lib/R/site-library:/usr/local/lib/R/library:/usr/lib/R/library'}" >> /usr/local/lib/R/etc/Renviron

## Set default CRAN repo
RUN echo 'options("repos"="http://cloud.r-project.org")' >> /usr/local/lib/R/etc/Rprofile.site

## More verbose UBSAN/SAN output (cf #3) -- this is still somewhat speculative
## Entry copied from Prof Ripley's setup described at http://www.stats.ox.ac.uk/pub/bdr/memtests/README.txt
ENV ASAN_OPTIONS 'alloc_dealloc_mismatch=0:detect_leaks=0:detect_odr_violation=0' 

RUN cd /usr/local/bin \
	&& mv R Rdevel \
	&& mv Rscript Rscriptdevel \
	&& ln -s Rdevel RD \
	&& ln -s Rscriptdevel RDscript



