* reorganize: move multistart explorations into a subdir?
* more investigation of bias?
* backporting to lme4pureR?
* explore remaining GH #936 issue for log-gaussian GLMMs: `misc/Gamma_GLMM/paramsurvey_loggaussian/` (B=100) confirms glmer's `sigma` for `gaussian(link="log")` is systematically biased low vs. glmmTMB/mgcv/glmmPQL (which agree tightly), not yet root-caused
* update the glmer vignette appropriately
