Package testing framework
===========================

Originally designed for `lme4`, but (??) should be extensible to other packages

Most of the additional jumping-through-hoops is intended to ensure that a complete local set of dependencies etc. etc. is installed ... might be better to set this up within a virtual machine.  Trying to get R only to look at a local library is a bit of a pain (see `setTestEnv`).

Next upgrade step would be to set up `make` rules so that packages were only tested when necessary; `Makefile` should require that the tarball be newer than the check `.Rout` files.

### To do

* current attempt 
 * failed to generate report because `R2HTML` wasn't available in the restricted `.libPaths`; 
 * tested against stable `lme4` because that was what was installed.  Fix.
 * didn't test some packages (`cplm`, `lmerTest`) because some dependencies were not installed/available (`biglm`, `ggplot2`, `reshape2`; `gplots`).  Why weren't these found?
* Write a makefile to test only tarballs newer than `check/PKG.Rcheck`
