\name{isNested}
\alias{isNested}
\title{Is f1 nested within f2?}
\usage{
  isNested(f1, f2)
}
\arguments{
  \item{f1}{factor 1}

  \item{f2}{factor 2}
}
\value{
  TRUE if factor 1 is nested within factor 2
}
\description{
  Does every level of f1 occur in conjunction with exactly
  one level of f2? The function is based on converting a
  triplet sparse matrix to a compressed column-oriented
  form in which the nesting can be quickly evaluated.
}
\examples{
with(Pastes, isNested(cask, batch))   ## => FALSE
with(Pastes, isNested(sample, batch))  ## => TRUE
}

