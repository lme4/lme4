* reorganize: move multistart explorations into a subdir?
* more investigation of bias?
* backporting to lme4pureR?
* check against GH #936, GH #557 -- GH #936 partially done: `misc/Gamma_GLMM/paramsurvey_loggaussian/` (B=100) confirms glmer's `sigma` for `gaussian(link="log")` is systematically biased low vs. glmmTMB/mgcv/glmmPQL (which agree tightly), not yet root-caused
* after merging master back into this branch, update the glmer vignette appropriately
