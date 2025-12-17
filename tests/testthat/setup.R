######################################################################
# Below is code frequently used for testing
# see: test-covariance_structures.R and test-covariance_nlmer.R

all.equal.nocheck <- function(x, y, ..., check.attributes = FALSE, check.class = FALSE) {
  ## working around mode-matching headaches
  if (is(x, "Matrix")) x <- matrix(x)
  if (is(y, "Matrix")) y <- matrix(y)
  all.equal(x, y, ..., check.attributes = check.attributes, check.class = check.class)
}

## set default tolerance to 5e-5 since we mostly use that
## 'tolerance' must be written out in full since it comes after ...
expect_equal_nocheck <- function(...,  tolerance = 5e-5) {
  expect_true(isTRUE(all.equal.nocheck(..., tolerance = tolerance)))
}

## Getting all equal as a number (in the all.equal examples documentation;
## don't know why they didn't make an argument instead!?)
all.eqNum <- function(...) {
  an <- all.equal.nocheck(...)
  if (isTRUE(an)) return(0)
  ## if check is less than tolerance all.equal returns TRUE, so sub() coerces to "TRUE"
  ##  and as.numeric() returns NA ...
  as.numeric(sub(".*:", '', an))
}
