Brain dump on interface issues
=====

How should we define the user interface for future/more flexible versions of `lme4`?

### All-in-one vs. little bits

What are the advantages of the all-in-one formula (i.e. `formula=y~x+(1|g)` rather than `formula=y~x, random=~1|g`)?

* It seems clearer (Doug likes it, and SW used to), but it seems to lead to additional code complexity for pulling apart/specifying e.g. whether you want the whole formula or just the fixed part or just the random part.  If we parsed the formula up-front (e.g. ran `findbars()/nobars()`) and stored the pieces separately, that would be a natural precursor to allowing flexibility in the interface (i.e, where users could use either interface style.
* The hard part in the specification is not necessarily whether fixed and random are separate but the specification of the random hard (which is intrinsically hard, as is the specification of the fixed part -- but people have more practice with that)
* Other than the fact the `nlme` and `lme4` interfaces are different, BMB doesn't feel like people have a particularly hard time with "oh, these are all supposed to be in the same formula" -- the hard part is figuring the distinction between an effect (i.e. LHS of the bar) and a grouping variable, the difference between single and multiple RE, the way that constructing a model matrix from the effect doesn't always work the way you think it should (e.g. the `dummy()` function)
* With just fixed + random specifications, it's not really a big deal whether we lump them together or not, but as we make the possibilities more complex, it gets to be a bigger and bigger deal

## Existing interfaces:

* SAS
* Stata (don't know much about this)?
* Genstat/AS-REML (not very sophisticated in GLMM world, still restricted to PQL), but very wide spectrum of linear model specification available due to the original ag experiment audience
* `MCMCglmm`: this may look very similar to the AS-REML interface. BMB finds MCMCglmm somewhat confusing, but maybe the problem lack of familiarity.
* A table of examples/correspondences would be nice.

## Interface issues for more complex models
* need to allow distance or coordinate metrics to be carried along with the data in cases where they're not reducible to a set of data frame columns -- simplest way is probably to attach attributes to the data frame and write accessor methods for attaching and extracting 
    * Could extend `data.frame` to allow additional things (probably as attributes? could this be done with a list?)
    * Alternatively could just make data a list with a one-line internal wrapper to embed data in a list if it's a data frame
   * transparency vs simplicity (list is most transparent, a new class is simplest)
* extra arguments: the value of `d(x|g,iid=TRUE)` vs. `d(x|g,iid=FALSE)` (equivalently just `d(x|g)` to specify a diagonal model with homogeneous vs. heterogeneous variance is debatable; `MCMCglmm` uses `idh()` vs `idv()` (I think) for the same distinction.  However, the additional arguments become essential for use cases like specifying variables and functions with which to compute distances.

## How much lme to use?

`lme` already has mature structures for specifying correlations, heteroscedasticity, and variance-covariance model structures.  How much of them should we (re-) use?

### Con:
* have to re-do some of Fabian's work
* many-separate-arguments structure 
* underlying design is baroque
* specification of distance computation is inflexible (e.g. it's not possible to make a Brownian motion model from the existing corExp structure because corExp insists on being x,y coordinates rather than allowing a raw distance matrix) -- assumes Euclidean distances
* current implementation is sparsely documented (we could fix this!)
* drops quickly into C code (not bad, but it would be nice to have a parallel R implementation available for inspection; corPhylo classes from ape are much easier to understand!  
* Perhaps too many levels of structure (again might be fixed by more documentation)?

### Pro:
* familiarity of user base
* several existing correlation structures built on lme (phylogenetic, great-circle, etc.)
* lots of fairly sophisticated code inside lme for computing Cholesky and inverse-Cholesky factors; we could use this code

* Perhaps, if we decide to use a different structure, we could write a wrapper that *converts* `corStruct`s to our preferred structure??

## Correlation structures

* idea of a *class* for correlation structures (as opposed to a function) is a good idea
* at least, it probably needs to contain separate functions for computing lambda and for processing theta into a meaningful format (currently done by embedding these functions in the correlation function's environment via `local()`, but I think this will be considerably harder for average users)
* take inspiration from `family()` -- allow a generic sort of constructor for correlation classes (but try not to make weird arbitrary-seeming design choices ...)
* SW suggested `corStruct.skeleton` to provide a starting point, but I think this will be achieved by an appropriate `corFamily()` (or whatever) function: with no arguments specified (i.e. all defaults), it would return something like a list of functions that would do the appropriate processing for a traditional unstructured (`(x|g)`) random effect. Structures like AR1, compound symmetry, would be special cases: roughly speaking, `ar1` would be to `corFamily()` as `binomial` is to `family()`. We might instead treat `corFamily()` as a constructor of class objects rather than functions, but the key point is that the objects thus constructed should be as *transparent* as possible!
* consider pdClasses as well -- how do these work, can we come up with a single class that combines the functionality of corStruct and pdClasses (and perhaps varStruct?)  In principle, the structures we want to be able to impose at the residual level and at the random effects level are more or less the same -- the key difference is whether we can use the block-diagonal structure of the residual variance (here we need to leverage semantic information about nested structures, or `isNested`, if possible!)



