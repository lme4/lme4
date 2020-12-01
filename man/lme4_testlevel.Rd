\name{lme4_testlevel}
\alias{lme4_testlevel}
\title{Detect testing level for lme4 examples and tests}
\description{
  Reads the environment variable \code{LME4_TEST_LEVEL} to determine
  which tests and examples to run
}
\usage{
lme4_testlevel()
}
\value{
  a numeric value: 1 for standard/'light' testing, larger values for
  more testing. Defaults to 1 if the environment variable is not set.
}

