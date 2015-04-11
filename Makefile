# from Yihui Xie, knitr package, but stripped down for now
# prepare the package for release
# NOTA BENE: This is *not* part of the CRAN release ==> can have GNUisms
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

build:
	cd ..;\
	 ## R CMD build --no-manual --no-build-vignettes $(PKGSRC)
	R CMD build $(PKGSRC)

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

travis: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --no-manual

integration-need:

	 ## git clone https://github.com/${TRAVIS_REPO_SLUG}-examples.git
	 ##	cd knitr-examples && \
	 ##	git checkout ${TRAVIS_BRANCH} && \
	 ##		GIT_PAGER=cat git show HEAD

integration-run: install
	rm knitr-examples/cache -rf
	make sysdeps deps xvfb-start knit xvfb-stop -C knitr-examples

integration-verify:
	GIT_PAGER=cat make diff -C knitr-examples

integration: integration-run integration-verify
