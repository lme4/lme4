
## references

- https://hackmd.io/W_7WkzZ_RRW249WHMlZ43Q
- https://bbolker.github.io/mixedmodels-misc/notes/varmats.html; https://github.com/bbolker/mixedmodels-misc/blob/master/notes/varmats.rmd
- https://github.com/rstats-gsoc/gsoc2025/wiki/Structured-covariance-matrices-for-lme4

There is nothing in the design of `lme4` that precludes us from implementing structured correlation matrices, by adding an additional step in the translation from $\theta$ parameters to $\Lambda$ (scaled Cholesky factor of the RE covariance matrix). We sketched out an approach years ago, but never brought it to completion.

* see also: https://github.com/lme4/lme4/issues/230, https://github.com/lme4/lme4/tree/flexLambda . We would probably have to sit down and reconstruct the design - I don't know if it's well-documented anywhere.

* strong opinions
  * `flexLambda` should be back-compatible with base `lme4` by design
  * write an update method for translating old structures to new, and/or include logic for accommodating old structures (or, write a new class? inheriting from old class?)
  
## examples/tests

* appropriate issues from `lme4` stack?
* implement some examples via modular interface (easy: diagonal, intermediate: CS (chol of CS is not trivial), AR1)
    * as wrapper
	* by hacking `Lambda`

## work

* start by going back to (close to) the branch point of `lme4` master/`flexLambda` and lift out the changes
	
## 

* allow alternative ways to construct `Lambda`
* hook for non-trivial functions to transform `theta` parameters to entries for `Lambda`

