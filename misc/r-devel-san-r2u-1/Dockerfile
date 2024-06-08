## Modified from https://gist.github.com/vjcitn/0cf0b0a25ec839ecbd8e0b47a23ce95e
FROM rocker/r-devel-san

RUN apt update -qq && apt install --yes --no-install-recommends ca-certificates gnupg
RUN wget -q -O- https://eddelbuettel.github.io/r2u/assets/dirk_eddelbuettel_key.asc     | tee -a /etc/apt/trusted.gpg.d/cranapt_key.asc
RUN echo "deb [arch=amd64] https://r2u.stat.illinois.edu/ubuntu jammy main"      > /etc/apt/sources.list.d/cranapt.list
RUN apt update -qq
RUN wget -q -O- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc     | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN echo "deb [arch=amd64] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"     > /etc/apt/sources.list.d/cran_r.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys     67C2D66C4B1D4339 51716619E084DAB9
RUN apt update -qq
## RUN ## Then install bspm (as root) and enable it, and enable a speed optimization
RUN install2.r bspm
RUN RHOME=$(RD RHOME)
RUN echo "suppressMessages(bspm::enable())" >> ${RHOME}/etc/Rprofile.site
RUN echo "options(bspm.version.check=FALSE)" >> ${RHOME}/etc/Rprofile.site
RUN echo "suppressMessages(bspm::enable())" >> ${HOME}/.Rprofile
RUN apt -y install r-cran-minqa r-cran-nloptr r-cran-rcpp r-cran-rcppeigen cmake
RUN echo "install.packages('lme4', dependencies = TRUE)" | RD --slave
RUN apt -y install git
RUN git clone https://github.com/eddelbuettel/sanitizers.git
RUN RD CMD INSTALL sanitizers
RUN wget https://cran.r-project.org/src/contrib/lme4_1.1-35.3.tar.gz
