Package testing framework
===========================

Originally designed for `lme4`, but (??) should be extensible to other packages

Most of the additional jumping-through-hoops is intended to ensure that a complete local set of dependencies etc. etc. is installed ... might be better to set this up within a virtual machine.  Trying to get R only to look at a local library is a bit of a pain (see `setTestEnv`).

HTML output of a relatively recent test can be viewed [here](http://htmlpreview.github.io/?https://github.com/lme4/lme4/blob/master/misc/pkgtests/lme4_compat_report.html)

### To do

* clean up "current version install" logic slightly?
* not sure dependencies-of-dependencies logic is correct ...
