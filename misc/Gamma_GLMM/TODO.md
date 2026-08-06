* reorganize: move multistart explorations into a subdir?
* more investigation of bias?
* backporting to lme4pureR?
* GH #936 log-gaussian `sigma` bias: mechanism now precisely characterized -- glmer's moment-based phi (`deviance/n` at the PIRLS conditional mode) is missing a `q` (number of random-effect levels) degrees-of-freedom correction; `phi_glmmTMB/phi_glmer` matches the predicted `n/(n-q)` to ~1% across B=100 replicates on the single-grouping-factor Rail example (see README). Not yet: (a) tried actually implementing an `n-q`-style correction to see if it closes the gap, (b) figured out what `q` should mean for crossed/multiple random-effect terms (glmmTMB handles these natively; `mgcv`/`glmmPQL`, the other two commensurate-`sigma` reference methods, don't -- limits how this could be validated beyond the single-factor case)
* update the glmer vignette appropriately
