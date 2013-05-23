`lme4` status
========================================================

Version 2013-05-08 17:03:31

## CURRENT (2013 May)

Almost on the list is "how can we release the damn thing?"

* as soon as feasible, update installation instructions, status, etc. [the lme4 r-forge project page](http://lme4.r-forge.r-project.org/) (put recent binary versions on the `lme4` repository, etc.)
* segfault/issues with R 3.0.+ *seem* to be resolved: can only reproduce segfaults on SCW's machine, with `debug` branch, with build from source (even after a clean Github checkout), with `example(glmer)` (`lmer` is OK).  BMB can't reproduce any more, on either Ubuntu or MacOS.
 * double-check with DB and MM that they can't reproduce any more either (*after* various DB fixes)
 * ??? (blue-sky) implement clean/virtual-machine installations with various permutations of `install_github`, git clone and install from source, install from tarball, various versions of upstream packages (CRAN vs R-forge), etc. (This is not going to happen any time soon!)
 * more realistically: check on BMB's virtual Windows machine??
 * otherwise, as far as we know, this is resolved -- probably want to put out some sort of vague warning / "let us know if you see anything fishy" to `r-sig-mixed-models`
 
* GLMM deviance function / PWRSS issues
 * further digging
 * ask DB whether state of deviance function and environments/reference classes associated therewith, are supposed to remain unchanged between successive calls.  If so, is this true only for a single optimization run (i.e. within the initial fitting process)?  It seems dangerous to give users a stateful deviance function (!!) BMB thinks (but may be wrong) that the purpose of maintaining state, other than the stuff that should be read-only (i.e. big data frames and matrices), would be to transmit information between successive steps of the PWRSS minimization process ... It might be inelegant, but should be easy, to reset values of things that should be clean at the beginning of a devfun call, either to zero or some sensible starting value, or to their values as they were initially set.  As long as we know that we're not breaking the design in this way, it should be (????) fairly easy to hack things right.  Could argue that we might be papering over cracks this way (rather than understanding what might be a problem) ...
  * collect any further test cases
  * as far as BMB recalls, brute-force resetting fixes every problem with possible exception of some very non-canonical links (which didn't work before anyway! and therefore are *not* release-critical)
  
* downstream packages: general
 * finish getting `pkgtests` working, follow up on issues with downstream developers (MM said long ago that he would be willing to help with this)

* downstream packages: `blme`
 * **very** short-term fix is to modify `blme` to look for `lme4.0` rather than `lme4` (i.e. search and replace `lme4.0` -> `lme4`)
 * lots of medium-term stuff that we hope we can do soon (see other notes ...)
 * highest priority: priors on variance-covariance matrices, and `glmer` with `family="gaussian"` to get (slow) priors on $\beta$
 
* **PAPERS**
 * PIRLS derivation; related to `glmer`/`blme` issues
 * do it!
 
## OLD (2013 Jan)

### Targets

* next two weeks (9-23 Jan): finish current writeups (GLMM procedure, RZX problem); finish BIRS proposal; maybe tackle another GLMM or NLMM issue (PWRSS failures?); maybe take a first shot at modularizing
* 16 Jan or earlier: draft of BIRS proposal to MM/DB
* 23 Jan or soon thereafter?: strategic Skype meeting with MM/DB

### Publishable units
* release of `lme4` unstable: ideally with working (non-fragile) GLMMs, but possibly as a fallback position with only LMMs fully supported
 * fix GLMMs????
 * get it working with downstream packages
  * modularize model-fitting process
  * coordinate with downstream package maintainers for more general issues (etc.)
  * change or augment `devComp` to include a native-R implementation?
* JRSS paper "#1": mostly written.  Update as necessary for any changes in computational structure or model interfaces. Update with more detailed information on GLMM fitting process? (Add SCW as co-author??)
* JRSS paper "#2": unwritten. More information about API, extensions, tricks and techniques, etc. (profiling, bootstrapping, calculating local curvatures/Wald variance-covariance matrices, post-hoc MCMC sampling) -- see `lme4-extras` vignette (still on r-forge, not copied yet to github)
* book:
 * update for changes in package
 * ?? add stuff from JRSS #1 and #2 (e.g. GLMM structure)
 * blue-sky chapter??
 
### glmer problem documentation
 
 Preparatory to dumping things in DB/MM's lap and whining to them to fix RcppEigen so that it can handle them ...
 
 * **fully** document the crabs/RZX problem.
  * we have already (1) instrumented `lme4` to show that it fails at the initial Cholesky decomposition of "B" (i.e. constructing `RZX`) (2) written out and checked native-R construction of all of the bits leading up to `B` (3) shown that `Matrix::Cholesky` gives us a sensible answer for the decomposition.  We haven't (4) instrumented `lme4.0` to show it does the same thing (5) written a `base::chol` example to show that it gets the same answer as `Matrix::Cholesky (6) run `RcppEigen::solveInPlace` in a minimal code file outside of `lme4` (7) do a thorough writeup of what's going on, including a `.RData` file that has all the pieces saved in it (i.e. original data, `Z`, `X`, `B`, `RZX`, ...) [#4 and #6 are **not necessary**]
* hopefully, ideally, if we can do it: diagnose and document an additional example (i.e. a PWRSS example)
 * maybe use an `nlmer` example?  (algal example)
 
### Level definitions

We found this useful in discussing `lme4` issues

* **Level I**: low-level computational/matrix representation stuff, i.e. what DB loves.  Julia, new `lme4` column-structured-branch stuff, etc.
* **Level II**: PWRSS/step-halving etc.
* **Level III**: optimization over $\theta$ and possibly $\beta$, using `nlminb`, `bobyqa`, etc. etc.. (e.g. Koller boundary-sticking problems)
* **Level IV**: API, formula interface, profiling, modularization of $Z$, etc.
